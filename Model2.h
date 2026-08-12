#pragma once

#include "Types.h"

#include <vector>
#include <unordered_map>
#include <functional>
#include <cstdint>
#include <utility>

class Model2 {
public:
    int _k;

    std::vector<InstanceID> _instances;
    NeibList _adj;

    const std::vector<FeatureType>& _features;

    std::unordered_map<FeatureType, int> FeatureNumber;
    std::vector<int> _featureIdOfInstance;

    std::vector<Color> ColorMap;

    std::vector<std::unordered_map<FeatureType, int>> InsDegree;

    std::vector<std::vector<InstanceID>> VSets;

    int size = 0;

    std::function<void(const std::vector<InstanceID>&)> onClique;
    std::vector<int> _edgeRank;
    std::vector<std::vector<InstanceID>> _kCliques;

public:
    Model2(const std::vector<InstanceID>& instances,
           const NeibList& insNeibs,
           const std::vector<FeatureType>& features,
           int k);

    void execute();

private:
    void buildFeatureStats();
    void buildColorMap();

    bool isValid(const std::vector<InstanceID>& curV, int need,
                 std::vector<int>& colorSeen,
                 std::vector<int>& featSeen,
                 int stamp) const;

    static void intersectVertices(const std::vector<InstanceID>& A,
                                  const std::vector<InstanceID>& B,
                                  std::vector<InstanceID>& out);

    void emitClique(const std::vector<InstanceID>& S);

    void collectActiveNeighbors(int eid,
                                const std::vector<std::vector<uint8_t>>& alive,
                                std::vector<InstanceID>& out) const;

private:
    struct SimpleEdge {
        InstanceID u, v;
    };

    std::vector<SimpleEdge> _edges;
    std::vector<std::vector<InstanceID>> _edgeVSets;

    // active masks：不物理删除，只标记当前是否有效
    std::vector<std::vector<uint8_t>> _edgeAlive;

    std::vector<int> _edgeOrder;
    std::unordered_map<std::uint64_t, int> _edgeId;
    std::vector<std::vector<int>> _incidentEdges;

    // 优化2：邻接边ID表，与_adj同步排序，用于快速获取边ID
    std::vector<std::vector<int>> _adjEdgeId;

    // 时间戳标记，避免递归热点中的 binary_search
    std::vector<int> _vertexMark;
    int _vertexMarkStamp = 0;

    static std::uint64_t edgeKey(InstanceID a, InstanceID b);

    void buildEdgeStructures();
    void buildIncidentEdges();
    void initEdgeAliveMasks();
    void EdgeDecomposition();

    bool markInactive(std::vector<std::vector<uint8_t>>& alive,
                      int eid,
                      InstanceID x) const;

    bool markInactiveWithLog(std::vector<std::vector<uint8_t>>& alive,
                             int eid,
                             InstanceID x,
                             std::vector<std::pair<int, int>>& changeLog) const;

    void dfsEBBkC(std::vector<InstanceID>& S,
                  std::vector<uint8_t>& usedFeat,
                  const std::vector<InstanceID>& cand,
                  const std::vector<int>& candESet,
                  int need,
                  int lastPos,
                  std::vector<int>& colorSeen,
                  std::vector<int>& featSeen,
                  int& stamp);
};