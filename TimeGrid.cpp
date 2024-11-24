#include "TimeGrid.h"
#include "analyticFunctions.h"

// constructor
TimeGrid::TimeGrid(
    const int _n1d, const double _dt, const double _final_time,
    const std::vector<std::vector<double>> &_terminal_value,
    const std::vector<std::vector<double>> &_obstacles) :
        n1d(_n1d), h(1/static_cast<double>(_n1d)), 
        dt(_dt), step_num(static_cast<int>(_final_time / _dt) + 1),
        grid(n1d, std::vector<std::vector<double>>(
             n1d, std::vector<double>(step_num, INFTY))) {
    // save obstacle boundaries
    for (const auto &obs : _obstacles) {
        std::vector<int> obs_index(4, 0);
        for (int i = 0; i < 4; ++i) {  
            obs_index[i] = static_cast<int>(obs[i] / h);
        }
        this->obstacles.push_back(obs_index);
    }
    // set terminal condition
    for (int i = 0; i < n1d; ++i) {
        for (int j = 0; j < n1d; ++j) {
            grid[i][j][step_num - 1] = _terminal_value[i][j];
        }
    }
}

// destructor
TimeGrid::~TimeGrid() { }

// Check if (x, y) is inside one of the obstacles
// Obstacles are considered as open sets
// Pay attention especially when obstacle 
// boundaries adhere to the domain boundaries!!!
bool TimeGrid::inObstacle(int j, int i) {
    for (const auto &obs_index : obstacles) {
        if (j > obs_index[0] && j < obs_index[2] &&
            i > obs_index[1] && i < obs_index[3]) {
            return true;
        }
    }
    return false;
}

// March backward in time starting from terminal condition
void TimeGrid::timeMarching() {
    for (int k = step_num - 1; k > 0; --k) {
        for (int i = 0; i < n1d; ++i) {
            for (int j = 0; j < n1d; ++j) {
                if (inObstacle(j, i))
                    continue;
                double x = static_cast<double>(j) * h;  
                double y = static_cast<double>(i) * h;

                double u_prev = grid[i][j][k];
                double f_prev = f(x, y);
                double K_prev = K(x, y);

                double u_x = upwindXGradient(i, j, k, u_prev);
                double u_y = upwindYGradient(i, j, k, u_prev);

                grid[i][j][k - 1] =  
                    u_prev - f_prev * sqrt(u_x * u_x + u_y * u_y) * dt +
                    K_prev * dt;
            }
        }
    }
}

double upwindScheme(double mid_right, double mid_left) {
    if (mid_right >= mid_left && mid_right > 0.0) {
        return mid_right;
    } else if (mid_left >= mid_right && mid_left > 0.0) {
        return mid_left;
    } else {
        return 0.0;
    }
}

double TimeGrid::upwindXGradient(int i, int j, int k, double u_prev) {
    double mid_right = (j == n1d - 1)
        ? 0.0
        : (u_prev - grid[i][j + 1][k]) / h;

    double mid_left = (j == 0)
        ? 0.0
        : (u_prev - grid[i][j - 1][k]) / h;
    return upwindScheme(mid_right, mid_left);
}
double TimeGrid::upwindYGradient(int i, int j, int k, double u_prev) {
    double mid_right = (i == n1d - 1)
        ? 0.0
        : (u_prev - grid[i + 1][j][k]) / h;

    double mid_left = (i == 0)
        ? 0.0
        : (u_prev - grid[i - 1][j][k]) / h;
    return upwindScheme(mid_right, mid_left);
}
