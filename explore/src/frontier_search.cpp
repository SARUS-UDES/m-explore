#include <explore/frontier_search.h>

#include <mutex>

#include <tf/tf.h>
#include <costmap_2d/cost_values.h>
#include <costmap_2d/costmap_2d.h>
#include <geometry_msgs/Point.h>

#include <explore/costmap_tools.h>

namespace frontier_exploration
{
using costmap_2d::LETHAL_OBSTACLE;
using costmap_2d::NO_INFORMATION;
using costmap_2d::FREE_SPACE;

class FrontierNormal {
public:
  enum direction : unsigned int {
    UP, DOWN, LEFT, RIGHT
  };

  FrontierNormal() {
    memset(normals, 0, sizeof(normals));
  }

  void accumulate(unsigned int x1, unsigned int y1, unsigned int x2, unsigned int y2) {
    FrontierNormal::direction dir;
    if (x1 < x2) {
      dir = FrontierNormal::direction::LEFT;
    }
    else if (x1 > x2) {
      dir = FrontierNormal::direction::RIGHT;
    }
    else if (y2 < y1) {
      dir = FrontierNormal::direction::UP;
    }
    else if (y2 > y1) {
      dir = FrontierNormal::direction::DOWN;
    }

    // Opposite directions cancel each other.
    // Normals cannot be both up and down or left and right simultaneously.
    const unsigned int opposite = direction_opposite[dir];
    if (normals[opposite] != 0) {
      normals[opposite]--;
    }
    else {
      normals[dir]++;
    }
  }

  double getAngle() const {
    const unsigned int nb_normals = normals[UP] + normals[DOWN] + normals[LEFT] + normals[RIGHT];
    double output = 0.0;
    if (nb_normals == 0) { return output; }

    // If there is a down component to the normal, set the right component to 2*pi.
    // Otherwise, the right component is equal to zero.
    if (normals[DOWN]) {
      output += normals[RIGHT] * 2 * M_PI;
    }
    output = normals[UP] * M_PI_2 + normals[DOWN] * -M_PI_2 + normals[LEFT] * M_PI;
    return output / nb_normals;
  }
  
private:
  unsigned int direction_opposite[4] = {DOWN, UP, RIGHT, LEFT};
  unsigned int normals[4];
};

FrontierSearch::FrontierSearch(costmap_2d::Costmap2D* costmap,
                               double orientation_scale,
                               double potential_scale, double gain_scale,
                               double min_frontier_size,
                               double max_frontier_size)
  : costmap_(costmap)
  , orientation_scale_(orientation_scale)
  , potential_scale_(potential_scale)
  , gain_scale_(gain_scale)
  , min_frontier_size_(min_frontier_size)
  , max_frontier_size_(max_frontier_size)
{
}

std::vector<Frontier> FrontierSearch::searchFrom(geometry_msgs::Pose pose)
{
  std::vector<Frontier> frontier_list;

  // Sanity check that robot is inside costmap bounds before searching
  unsigned int mx, my;
  if (!costmap_->worldToMap(pose.position.x, pose.position.y, mx, my)) {
    ROS_ERROR("Robot out of costmap bounds, cannot search for frontiers");
    return frontier_list;
  }

  double yaw;
  {
    tf::Quaternion q(
      pose.orientation.x,
      pose.orientation.y,
      pose.orientation.z,
      pose.orientation.w);
    tf::Matrix3x3 m(q);
    double roll, pitch;
    m.getRPY(roll, pitch, yaw);
  }

  // make sure map is consistent and locked for duration of search
  std::lock_guard<costmap_2d::Costmap2D::mutex_t> lock(*(costmap_->getMutex()));

  map_ = costmap_->getCharMap();
  size_x_ = costmap_->getSizeInCellsX();
  size_y_ = costmap_->getSizeInCellsY();

  // initialize flag arrays to keep track of visited and frontier cells
  std::vector<bool> frontier_flag(size_x_ * size_y_, false);
  std::vector<bool> visited_flag(size_x_ * size_y_, false);

  // initialize breadth first search
  std::queue<unsigned int> bfs;

  // find closest clear cell to start search
  unsigned int clear, pos = costmap_->getIndex(mx, my);
  if (nearestCell(clear, pos, FREE_SPACE, *costmap_)) {
    bfs.push(clear);
  } else {
    bfs.push(pos);
    ROS_WARN("Could not find nearby clear cell to start search");
  }
  visited_flag[bfs.front()] = true;

  while (!bfs.empty()) {
    unsigned int idx = bfs.front();
    bfs.pop();

    // iterate over 4-connected neighbourhood
    for (unsigned nbr : nhood4(idx, *costmap_)) {
      // add to queue all free, unvisited cells, use descending search in case
      // initialized on non-free cell
      if (map_[nbr] <= map_[idx] && !visited_flag[nbr]) {
        visited_flag[nbr] = true;
        bfs.push(nbr);
        // check if cell is new frontier cell (unvisited, NO_INFORMATION, free
        // neighbour)
      } else if (isNewFrontierCell(nbr, frontier_flag)) {
        frontier_flag[nbr] = true;
        Frontier new_frontier = buildNewFrontier(nbr, pos, frontier_flag);
        if (new_frontier.size * costmap_->getResolution() >=
            min_frontier_size_) {
          frontier_list.push_back(new_frontier);
        }
      }
    }
  }

  // set costs of frontiers
  for (auto& frontier : frontier_list) {
    frontier.cost = frontierCost(frontier, pose, yaw);
  }
  std::sort(
      frontier_list.begin(), frontier_list.end(),
      [](const Frontier& f1, const Frontier& f2) { return f1.cost < f2.cost; });

  return frontier_list;
}

Frontier FrontierSearch::buildNewFrontier(unsigned int initial_cell,
                                          unsigned int reference,
                                          std::vector<bool>& frontier_flag)
{
  const double costmap_res = costmap_->getResolution();
  const double max_size = max_frontier_size_ / costmap_res;

  // initialize frontier structure
  Frontier output;
  output.centroid.x = 0;
  output.centroid.y = 0;
  output.size = 1;
  output.min_distance = std::numeric_limits<double>::infinity();

  // record initial contact point for frontier
  unsigned int ix, iy;
  costmap_->indexToCells(initial_cell, ix, iy);
  costmap_->mapToWorld(ix, iy, output.initial.x, output.initial.y);

  // push initial gridcell onto queue
  std::queue<unsigned int> bfs;
  bfs.push(initial_cell);

  // cache reference position in world coords
  unsigned int rx, ry;
  double reference_x, reference_y;
  costmap_->indexToCells(reference, rx, ry);
  costmap_->mapToWorld(rx, ry, reference_x, reference_y);

  while (!bfs.empty()) {
    unsigned int idx = bfs.front();
    bfs.pop();

    // try adding cells in 8-connected neighborhood to frontier
    for (unsigned int nbr : nhood8(idx, *costmap_)) {
      // check if neighbour is a potential frontier cell
      if (isNewFrontierCell(nbr, frontier_flag)) {
        // mark cell as frontier
        frontier_flag[nbr] = true;
        unsigned int mx, my;
        double wx, wy;
        costmap_->indexToCells(nbr, mx, my);
        costmap_->mapToWorld(mx, my, wx, wy);

        geometry_msgs::Point point;
        point.x = wx;
        point.y = wy;
        output.points.push_back(point);

        // update frontier size
        output.size++;

        // update centroid of frontier
        output.centroid.x += wx;
        output.centroid.y += wy;

        // add to queue for breadth first search
        bfs.push(nbr);
        if (output.size >= max_size) {
          break;
        }
      }
    }
    if (output.size >= max_size) {
      break;
    }
  }

  // average out frontier centroid
  output.centroid.x /= output.size;
  output.centroid.y /= output.size;
  // Manhattan distance is good enough since this is used to rank frontiers based on distances.
  output.min_distance = abs(static_cast<double>(reference_x - output.centroid.x)) + abs(static_cast<double>(reference_y - output.centroid.y));
  output.middle.x = output.centroid.x;
  output.middle.y = output.centroid.y;
  return output;
}

double FrontierSearch::computeNormal(const std::vector<geometry_msgs::Point>& points)
{
  FrontierNormal frontier_normal;

  for (geometry_msgs::Point pt : points) {
    unsigned int mx, my;
    costmap_->worldToMap(pt.x, pt.y, mx, my);
    unsigned int idx = costmap_->getIndex(mx, my);
    // Compute normals
    for (unsigned int nbr: nhood4(idx, *costmap_)) {
      if (map_[nbr] == FREE_SPACE){
        unsigned int x1, y1, x2, y2;
        costmap_->indexToCells(idx, x1, y1);
        costmap_->indexToCells(nbr, x2, y2);
        frontier_normal.accumulate(x1, y1, x2, y2);
      }
    }
  }
  return frontier_normal.getAngle();;
}

bool FrontierSearch::isNewFrontierCell(unsigned int idx,
                                       const std::vector<bool>& frontier_flag)
{
  // check that cell is unknown and not already marked as frontier
  if (map_[idx] != NO_INFORMATION || frontier_flag[idx]) {
    return false;
  }

  // frontier cells should have at least one cell in 4-connected neighbourhood
  // that is free
  for (unsigned int nbr : nhood4(idx, *costmap_)) {
    if (map_[nbr] == FREE_SPACE) {
      return true;
    }
  }

  return false;
}

double FrontierSearch::frontierCost(const Frontier& frontier, const geometry_msgs::Pose& pose, double yaw)
{
  const double v2_x = pose.position.x - frontier.centroid.x;
  const double v2_y = pose.position.y - frontier.centroid.y;
  const double heading_with_frontier = (M_PI - abs(atan2f64(v2_y, v2_x) - yaw)) / M_PI;
  const double proximity_cost = costmap_->getResolution()*(potential_scale_ * frontier.min_distance - gain_scale_ * frontier.size/max_frontier_size_);
  return orientation_scale_*(heading_with_frontier) + proximity_cost;
}
}
