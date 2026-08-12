#pragma once

#include <vector>
#include <map>
#include "Types.h"

/**
 * 基于平面坐标的邻居查找：
 * - 输入：实例 ID 列表 + 各实例的坐标 & 特征
 * - 输出：InsNeibs 邻接表（InstanceID -> 邻居实例 ID 列表）
 *
 * 要求：
 *   1) 两个实例距离 <= _distance
 *   2) 两个实例的 feature 不相同
 */
class FindNeibs {
public:
    // 邻接表：InsNeibs[u] = 所有与 u 相邻的实例 ID
    NeibList InsNeibs;

private:
    // 网格 -> 该网格中的实例 ID 列表
    std::map<CellId, std::vector<InstanceID>> gridIns;

    // 所有实例的坐标 & 特征（由外部传入，并在类中仅持有引用，不拷贝）
    const std::vector<LocationType>& _locations;
    const std::vector<FeatureType>& _features;

    double _distance;
    int grid_X = 0;
    int grid_Y = 0;

public:
    /**
     * @param instanceIds  需要构图的实例 ID 列表（假定为连续 int，0..N-1 或子集）
     * @param locations    全局实例坐标数组：locations[id] = (x,y)
     * @param features     全局实例特征数组：features[id] = featureName
     * @param distance     判定邻居的最大欧氏距离
     * @param min_x        所有实例中 x 的最小值（用于网格编号）
     * @param min_y        所有实例中 y 的最小值
     * @param max_x        所有实例中 x 的最大值（目前仅用于统计最大网格编号，可选）
     * @param max_y        所有实例中 y 的最大值
     */
    FindNeibs(const std::vector<InstanceID>& instanceIds,
              const std::vector<LocationType>& locations,
              const std::vector<FeatureType>& features,
              double distance,
              double min_x, double min_y,
              double max_x, double max_y);

private:
    // 判断两个实例是否在距离阈值内
    bool isNeibs(InstanceID a, InstanceID b) const;

    // 获取某个网格的“自身 + 周围 8 个邻接格子”（只返回存在于 gridIns 的格子）
    std::vector<CellId> getNeibGrids(const CellId& centerGrid) const;

    // 核心：在网格结构上遍历所有实例对，填充 InsNeibs
    void getNeibs(const std::vector<InstanceID>& instanceIds);
};
