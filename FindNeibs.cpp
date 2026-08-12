#include "FindNeibs.h"
#include <cmath>
#include <algorithm>

FindNeibs::FindNeibs(const std::vector<InstanceID>& instanceIds,
                     const std::vector<LocationType>& locations,
                     const std::vector<FeatureType>& features,
                     double distance,
                     double min_x, double min_y,
                     double max_x, double max_y)
        : _locations(locations),
          _features(features),
          _distance(distance),
          grid_X(0),
          grid_Y(0)
{
    // 1) 根据坐标把实例丢进网格
    int maxId = 0;
    for (auto id : instanceIds) {
        maxId = std::max(maxId, id);
        const auto& loc = _locations[id];
        int grid_x = static_cast<int>((loc.first  - min_x) / distance);
        int grid_y = static_cast<int>((loc.second - min_y) / distance);
        CellId cell{grid_x, grid_y};
        gridIns[cell].push_back(id);

        if (grid_x > grid_X) grid_X = grid_x;
        if (grid_y > grid_Y) grid_Y = grid_y;
    }

    // 2) 初始化邻接表大小（按最大 ID 分配）
    InsNeibs.clear();
    InsNeibs.resize(maxId + 1);

    // 3) 在网格结构上计算邻居
    getNeibs(instanceIds);
}

bool FindNeibs::isNeibs(InstanceID a, InstanceID b) const {
    const auto& loc1 = _locations[a];
    const auto& loc2 = _locations[b];
    double dx = loc1.first  - loc2.first;
    double dy = loc1.second - loc2.second;
    double d2 = dx * dx + dy * dy;
    return d2 <= _distance * _distance;
}

std::vector<CellId> FindNeibs::getNeibGrids(const CellId& centerGrid) const {
    std::vector<CellId> grids;
    int gx = centerGrid.first;
    int gy = centerGrid.second;

    auto try_add = [&](int x, int y) {
        CellId c{x, y};
        if (gridIns.count(c)) grids.push_back(c);
    };

    // 自身 + 八邻格
    try_add(gx,     gy);
    try_add(gx - 1, gy - 1);
    try_add(gx,     gy - 1);
    try_add(gx + 1, gy - 1);
    try_add(gx - 1, gy);
    try_add(gx + 1, gy);
    try_add(gx - 1, gy + 1);
    try_add(gx,     gy + 1);
    try_add(gx + 1, gy + 1);

    return grids;
}

void FindNeibs::getNeibs(const std::vector<InstanceID>& instanceIds) {
    // 对每个网格，取自身 + 邻接网格中的实例，判定是否为邻居
    for (auto it = gridIns.begin(); it != gridIns.end(); ++it) {
        CellId centerGrid = it->first;
        const auto& insInCenter = it->second;

        std::vector<CellId> neibGrids = getNeibGrids(centerGrid);

        for (auto ins1 : insInCenter) {
            const auto& f1 = _features[ins1];

            for (const auto& cell : neibGrids) {
                const auto& insInCell = gridIns.at(cell);
                for (auto ins2 : insInCell) {
                    if (ins1 == ins2) continue;                  // 同一个实例跳过
                    if (f1 == _features[ins2]) continue;         // 特征相同跳过

                    // 为避免重复添加边，只在 ins1 < ins2 时建立一次
                    if (ins1 < ins2 && isNeibs(ins1, ins2)) {
                        InsNeibs[ins1].push_back(ins2);
                        InsNeibs[ins2].push_back(ins1);
                    }
                }
            }
        }
    }
}
