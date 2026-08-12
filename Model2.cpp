#include "Model2.h"

#include <algorithm>
#include <queue>
#include <vector>
#include <cstdint>
#include <unordered_map>

Model2::Model2(const std::vector<InstanceID>& instances,
               const NeibList& insNeibs,
               const std::vector<FeatureType>& features,
               int k)
        : _k(k),
          _instances(instances),
          _adj(insNeibs),
          _features(features)
{
    const int N = static_cast<int>(_features.size());
    ColorMap.assign(N, -1);
    InsDegree.assign(N, {});
    VSets.assign(N, {});
    _featureIdOfInstance.assign(N, -1);
    _vertexMark.assign(N, 0);
    buildFeatureStats();
}

void Model2::buildFeatureStats() {
    FeatureNumber.clear();
    FeatureNumber.reserve(_instances.size());

    int fid = 0;
    for (auto id : _instances) {
        const auto& f = _features[id];
        if (!FeatureNumber.count(f)) {
            FeatureNumber[f] = fid++;
        }
        _featureIdOfInstance[id] = FeatureNumber[f];
    }

    for (size_t u = 0; u < _adj.size(); ++u) {
        for (auto v : _adj[u]) {
            const auto& fv = _features[v];
            const auto& fu = _features[u];
            if (fu == fv) continue;

            ++InsDegree[u][fv];
            ++InsDegree[v][fu];
        }
    }
}

void Model2::buildColorMap() {
    std::vector<int> usedStamp(_features.size() + 1, 0);
    int stamp = 0;

    for (auto id : _instances) {
        ++stamp;
        for (auto nb : _adj[id]) {
            int c = ColorMap[nb];
            if (c >= 0) usedStamp[c] = stamp;
        }

        int c = 0;
        while (usedStamp[c] == stamp) ++c;
        ColorMap[id] = c;
    }
}

void Model2::buildEdgeStructures() {
    const int N = static_cast<int>(_features.size());
    std::vector<char> inCore(N, 0);
    for (auto id : _instances) inCore[id] = 1;

    _edges.clear();
    _edgeVSets.clear();
    _edgeOrder.clear();
    _edgeId.clear();

    for (auto u : _instances) {
        auto& list = _adj[u];
        std::sort(list.begin(), list.end());
        list.erase(std::unique(list.begin(), list.end()), list.end());
    }

    _edgeId.reserve(_instances.size() * 4);
    for (auto u : _instances) {
        for (auto v : _adj[u]) {
            if (!inCore[v]) continue;
            if (u < v) {
                int eid = static_cast<int>(_edges.size());
                _edges.push_back(SimpleEdge{u, v});
                _edgeVSets.emplace_back();
                _edgeId[edgeKey(u, v)] = eid;
            }
        }
    }

    for (int eid = 0; eid < static_cast<int>(_edges.size()); ++eid) {
        InstanceID u = _edges[eid].u;
        InstanceID v = _edges[eid].v;
        const auto& Nu = _adj[u];
        const auto& Nv = _adj[v];
        std::vector<InstanceID> tmp;
        intersectVertices(Nu, Nv, tmp);
        _edgeVSets[eid] = std::move(tmp);
    }

    // 构建 _adj 和同步的 _adjEdgeId
    _adj.assign(_features.size(), {});
    _adjEdgeId.assign(_features.size(), {});
    for (int eid = 0; eid < static_cast<int>(_edges.size()); ++eid) {
        const auto& e = _edges[eid];
        _adj[e.u].push_back(e.v);
        _adj[e.v].push_back(e.u);
        _adjEdgeId[e.u].push_back(eid);
        _adjEdgeId[e.v].push_back(eid);
    }
    // 排序去重并同步 _adjEdgeId
    for (size_t i = 0; i < _adj.size(); ++i) {
        auto& adjList = _adj[i];
        auto& idList = _adjEdgeId[i];
        if (adjList.empty()) continue;
        std::vector<std::pair<InstanceID, int>> paired(adjList.size());
        for (size_t j = 0; j < adjList.size(); ++j) {
            paired[j] = {adjList[j], idList[j]};
        }
        std::sort(paired.begin(), paired.end());
        adjList.clear();
        idList.clear();
        for (auto& p : paired) {
            if (adjList.empty() || adjList.back() != p.first) {
                adjList.push_back(p.first);
                idList.push_back(p.second);
            }
        }
    }
}

void Model2::buildIncidentEdges() {
    _incidentEdges.assign(_features.size(), {});
    for (int eid = 0; eid < (int)_edges.size(); ++eid) {
        const auto& e = _edges[eid];
        _incidentEdges[e.u].push_back(eid);
        _incidentEdges[e.v].push_back(eid);
    }
}

void Model2::initEdgeAliveMasks() {
    _edgeAlive.clear();
    _edgeAlive.resize(_edgeVSets.size());
    for (int eid = 0; eid < (int)_edgeVSets.size(); ++eid) {
        _edgeAlive[eid].assign(_edgeVSets[eid].size(), 1);
    }
}

// 辅助：两个有序向量的交集（输出参数版本）
void Model2::intersectVertices(const std::vector<InstanceID>& A,
                               const std::vector<InstanceID>& B,
                               std::vector<InstanceID>& out) {
    out.clear();
    out.reserve(std::min(A.size(), B.size()));
    size_t i = 0, j = 0;
    while (i < A.size() && j < B.size()) {
        if (A[i] < B[j]) ++i;
        else if (B[j] < A[i]) ++j;
        else {
            out.push_back(A[i]);
            ++i; ++j;
        }
    }
}

void Model2::EdgeDecomposition() {
    const int k = _k;
    const int need = k - 2;
    const int m = static_cast<int>(_edges.size());
    if (m == 0 || need <= 0) return;

    std::vector<char> edgeAlive(m, 1);

    // 优化1：使用 unordered_map 替代 vector<pair<int,int>>
    std::vector<std::unordered_map<int, int>> edgeColorCnt(m);
    std::vector<std::unordered_map<int, int>> edgeFeatCnt(m);
    std::vector<int> edgeColorKinds(m, 0);
    std::vector<int> edgeFeatKinds(m, 0);

    std::vector<std::vector<uint8_t>> cnAlive(m);

    // 优化2：使用 _adjEdgeId 二分查找代替 _edgeId 哈希表
    auto getEid = [&](InstanceID a, InstanceID b) -> int {
        const auto& adjList = _adj[a];
        auto it = std::lower_bound(adjList.begin(), adjList.end(), b);
        if (it != adjList.end() && *it == b) {
            int pos = it - adjList.begin();
            return _adjEdgeId[a][pos];
        }
        return -1;
    };

    auto isBadEdge = [&](int eid) -> bool {
        return edgeColorKinds[eid] < need || edgeFeatKinds[eid] < need;
    };

    // 在 unordered_map 中增加计数，并更新种类数
    auto incCount = [](std::unordered_map<int, int>& mp, int key, int& kinds) {
        auto it = mp.find(key);
        if (it == mp.end()) {
            mp[key] = 1;
            ++kinds;
        } else {
            ++(it->second);
        }
    };
    // 在 unordered_map 中减少计数，若减到0则删除，并更新种类数
    auto decCount = [](std::unordered_map<int, int>& mp, int key, int& kinds) -> bool {
        auto it = mp.find(key);
        if (it == mp.end()) return false;
        if (--(it->second) == 0) {
            mp.erase(it);
            --kinds;
            return true;
        }
        return false;
    };

    for (int eid = 0; eid < m; ++eid) {
        const auto &cn = _edgeVSets[eid];
        cnAlive[eid].assign(cn.size(), 1);

        auto &cCnt = edgeColorCnt[eid];
        auto &fCnt = edgeFeatCnt[eid];
        // 预分配空间（可选，unordered_map 不支持 reserve 精确值，但可以给个大致值）
        cCnt.reserve(std::min(cn.size(), size_t(256)));
        fCnt.reserve(std::min(cn.size(), size_t(256)));

        int cKinds = 0, fKinds = 0;
        for (auto w : cn) {
            incCount(cCnt, (int)ColorMap[w], cKinds);
            incCount(fCnt, _featureIdOfInstance[w], fKinds);
        }
        edgeColorKinds[eid] = cKinds;
        edgeFeatKinds[eid]  = fKinds;
    }

    auto removeCommonNeighbor = [&](int eid, InstanceID x) -> bool {
        auto &cn = _edgeVSets[eid];
        auto &al = cnAlive[eid];

        auto it = std::lower_bound(cn.begin(), cn.end(), x);
        if (it == cn.end() || *it != x) return false;

        const int pos = (int)(it - cn.begin());
        if (al[pos] == 0) return false;
        al[pos] = 0;

        int cx = (int)ColorMap[x];
        decCount(edgeColorCnt[eid], cx, edgeColorKinds[eid]);

        int fx = _featureIdOfInstance[x];
        decCount(edgeFeatCnt[eid], fx, edgeFeatKinds[eid]);

        return true;
    };

    std::queue<int> Q;
    std::vector<char> inQ(m, 0);

    auto pushIfNeeded = [&](int eid) {
        if (eid < 0 || eid >= m) return;
        if (!edgeAlive[eid]) return;
        if (!inQ[eid]) {
            Q.push(eid);
            inQ[eid] = 1;
        }
    };

    for (int eid = 0; eid < m; ++eid) {
        if (isBadEdge(eid)) pushIfNeeded(eid);
    }

    while (!Q.empty()) {
        int eid = Q.front();
        Q.pop();
        inQ[eid] = 0;

        if (!edgeAlive[eid]) continue;
        if (!isBadEdge(eid)) continue;

        edgeAlive[eid] = 0;
        const InstanceID u = _edges[eid].u;
        const InstanceID v = _edges[eid].v;

        const auto &cn = _edgeVSets[eid];
        auto &al = cnAlive[eid];

        for (int i = 0; i < (int)cn.size(); ++i) {
            if (!al[i]) continue;
            InstanceID w = cn[i];

            int e_uw = getEid(u, w);
            if (e_uw != -1 && edgeAlive[e_uw]) {
                if (removeCommonNeighbor(e_uw, v) && isBadEdge(e_uw)) {
                    pushIfNeeded(e_uw);
                }
            }

            int e_vw = getEid(v, w);
            if (e_vw != -1 && edgeAlive[e_vw]) {
                if (removeCommonNeighbor(e_vw, u) && isBadEdge(e_vw)) {
                    pushIfNeeded(e_vw);
                }
            }
        }

        cnAlive[eid].clear();
        edgeColorCnt[eid].clear();
        edgeFeatCnt[eid].clear();
        edgeColorKinds[eid] = 0;
        edgeFeatKinds[eid]  = 0;
    }

    std::vector<SimpleEdge> newEdges;
    std::vector<std::vector<InstanceID>> newEdgeVSets;
    newEdges.reserve(m);
    newEdgeVSets.reserve(m);

    for (int eid = 0; eid < m; ++eid) {
        if (!edgeAlive[eid]) continue;

        newEdges.push_back(_edges[eid]);

        auto &cn = _edgeVSets[eid];
        auto &al = cnAlive[eid];

        std::vector<InstanceID> compact;
        compact.reserve(cn.size());
        for (int i = 0; i < (int)cn.size(); ++i) {
            if (al[i]) compact.push_back(cn[i]);
        }
        newEdgeVSets.push_back(std::move(compact));
    }

    _edges.swap(newEdges);
    _edgeVSets.swap(newEdgeVSets);

    _edgeId.clear();
    _edgeId.reserve(_edges.size() * 2 + 1);
    for (int eid = 0; eid < (int)_edges.size(); ++eid) {
        _edgeId[edgeKey(_edges[eid].u, _edges[eid].v)] = eid;
    }

    _edgeOrder.resize(_edges.size());
    for (int i = 0; i < (int)_edges.size(); ++i) _edgeOrder[i] = i;

    std::sort(_edgeOrder.begin(), _edgeOrder.end(),
              [&](int a, int b) { return _edgeVSets[a].size() < _edgeVSets[b].size(); });

    _edgeRank.assign(_edges.size(), 0);
    for (int pos = 0; pos < (int)_edgeOrder.size(); ++pos) {
        _edgeRank[_edgeOrder[pos]] = pos;
    }

    // 重建 _adj 和同步的 _adjEdgeId
    _adj.assign(_features.size(), {});
    _adjEdgeId.assign(_features.size(), {});
    for (int eid = 0; eid < (int)_edges.size(); ++eid) {
        const auto& e = _edges[eid];
        _adj[e.u].push_back(e.v);
        _adj[e.v].push_back(e.u);
        _adjEdgeId[e.u].push_back(eid);
        _adjEdgeId[e.v].push_back(eid);
    }
    // 排序去重并同步 _adjEdgeId
    for (size_t i = 0; i < _adj.size(); ++i) {
        auto& adjList = _adj[i];
        auto& idList = _adjEdgeId[i];
        if (adjList.empty()) continue;
        std::vector<std::pair<InstanceID, int>> paired(adjList.size());
        for (size_t j = 0; j < adjList.size(); ++j) {
            paired[j] = {adjList[j], idList[j]};
        }
        std::sort(paired.begin(), paired.end());
        adjList.clear();
        idList.clear();
        for (auto& p : paired) {
            if (adjList.empty() || adjList.back() != p.first) {
                adjList.push_back(p.first);
                idList.push_back(p.second);
            }
        }
    }
}

bool Model2::isValid(const std::vector<InstanceID>& curV, int need,
                     std::vector<int>& colorSeen,
                     std::vector<int>& featSeen,
                     int stamp) const {
    if ((int)curV.size() < need) return false;

    int colorKinds = 0;
    int featKinds = 0;
    for (auto id : curV) {
        const int c = ColorMap[id];
        if (colorSeen[c] != stamp) {
            colorSeen[c] = stamp;
            ++colorKinds;
        }

        const int f = _featureIdOfInstance[id];
        if (featSeen[f] != stamp) {
            featSeen[f] = stamp;
            ++featKinds;
        }

        if (colorKinds >= need && featKinds >= need) return true;
    }

    return colorKinds >= need && featKinds >= need;
}

void Model2::emitClique(const std::vector<InstanceID>& S) {
    if (onClique) {
        onClique(S);
    }
}

// 输出参数版本：收集活跃公共邻居
void Model2::collectActiveNeighbors(int eid,
                                    const std::vector<std::vector<uint8_t>>& alive,
                                    std::vector<InstanceID>& out) const {
    out.clear();
    const auto& base = _edgeVSets[eid];
    const auto& mask = alive[eid];
    out.reserve(base.size());
    for (int i = 0; i < (int)base.size(); ++i) {
        if (mask[i]) out.push_back(base[i]);
    }
}

// 使用二分查找代替 _edgePos
bool Model2::markInactive(std::vector<std::vector<uint8_t>>& alive,
                          int eid,
                          InstanceID x) const {
    if (eid < 0 || eid >= (int)alive.size()) return false;
    const auto& base = _edgeVSets[eid];
    auto it = std::lower_bound(base.begin(), base.end(), x);
    if (it == base.end() || *it != x) return false;
    const int pos = (int)(it - base.begin());
    if (alive[eid][pos] == 0) return false;
    alive[eid][pos] = 0;
    return true;
}

bool Model2::markInactiveWithLog(std::vector<std::vector<uint8_t>>& alive,
                                 int eid,
                                 InstanceID x,
                                 std::vector<std::pair<int,int>>& changeLog) const {
    if (eid < 0 || eid >= (int)alive.size()) return false;
    const auto& base = _edgeVSets[eid];
    auto it = std::lower_bound(base.begin(), base.end(), x);
    if (it == base.end() || *it != x) return false;
    const int pos = (int)(it - base.begin());
    if (alive[eid][pos] == 0) return false;
    alive[eid][pos] = 0;
    changeLog.emplace_back(eid, pos);
    return true;
}

void Model2::dfsEBBkC(std::vector<InstanceID>& S,
                      std::vector<uint8_t>& usedFeat,
                      const std::vector<InstanceID>& cand,
                      const std::vector<int>& candESet,
                      int need, int lastPos,
                      std::vector<int>& colorSeen,
                      std::vector<int>& featSeen,
                      int& stamp)
{
    if ((int)cand.size() < need || (int)candESet.size() < need*(need-1)/2) return;

    ++stamp;
    if (!isValid(cand, need, colorSeen, featSeen, stamp)) return;

    if (need >= 2 && candESet.empty()) return;
    if(candESet.size()< need*(need-1)/2) return;

    if (need == 0) {
        emitClique(S);
        size++;
        return;
    }

    if (need == 1) {
        for (auto x : cand) {
            int fx = _featureIdOfInstance[x];
            if (usedFeat[fx]) continue;

            S.push_back(x);
            usedFeat[fx] = 1;

            emitClique(S);
            size++;

            usedFeat[fx] = 0;
            S.pop_back();
        }
        return;
    }

    if (need == 2) {
        const int m = (int)candESet.size();
        int start = 0;
        {
            int lo = 0, hi = m;
            while (lo < hi) {
                int mid = (lo + hi) >> 1;
                if (_edgeRank[candESet[mid]] <= lastPos) lo = mid + 1;
                else hi = mid;
            }
            start = lo;
        }

        for (int i = start; i < m; ++i) {
            const int eid = candESet[i];
            const InstanceID u = _edges[eid].u;
            const InstanceID v = _edges[eid].v;

            const int fu = _featureIdOfInstance[u];
            const int fv = _featureIdOfInstance[v];

            if (fu == fv) continue;
            if (usedFeat[fu] || usedFeat[fv]) continue;

            S.push_back(u);
            S.push_back(v);

            emitClique(S);
            size++;

            S.pop_back();
            S.pop_back();
        }
        return;
    }

    const int m = (int)candESet.size();
    int start = 0;
    {
        int lo = 0, hi = m;
        while (lo < hi) {
            int mid = (lo + hi) >> 1;
            if (_edgeRank[candESet[mid]] <= lastPos) lo = mid + 1;
            else hi = mid;
        }
        start = lo;
    }

    std::vector<std::pair<int,int>> frameChangeLog;
    frameChangeLog.reserve(256);

    // 复用缓冲区
    std::vector<InstanceID> activeNbrs;
    std::vector<InstanceID> nextCand;

    for (int i = start; i < m; ++i) {
        const int eid = candESet[i];
        const int newLastPos = _edgeRank[eid];

        const InstanceID u = _edges[eid].u;
        const InstanceID v = _edges[eid].v;

        const int fu = _featureIdOfInstance[u];
        const int fv = _featureIdOfInstance[v];

        if (fu == fv) continue;
        if (usedFeat[fu] || usedFeat[fv]) continue;

        S.push_back(u);
        S.push_back(v);
        usedFeat[fu] = 1;
        usedFeat[fv] = 1;

        // 收集活跃公共邻居（复用缓冲区）
        collectActiveNeighbors(eid, _edgeAlive, activeNbrs);

        // Nextcand = cand ∩ activeNbrs
        intersectVertices(activeNbrs, cand, nextCand);

        if ((int)nextCand.size() >= need - 2) {
            std::vector<int> nextCandESet;
            nextCandESet.reserve(m - (i + 1));

            ++_vertexMarkStamp;
            const int mark = _vertexMarkStamp;
            for (InstanceID x : nextCand) _vertexMark[x] = mark;

            for (int j = i + 1; j < m; ++j) {
                const int e2 = candESet[j];
                const InstanceID a = _edges[e2].u;
                const InstanceID b = _edges[e2].v;

                if (_vertexMark[a] != mark) continue;
                if (_vertexMark[b] != mark) continue;
                if (_featureIdOfInstance[a] == _featureIdOfInstance[b]) continue;

                nextCandESet.push_back(e2);
            }

            dfsEBBkC(S, usedFeat, nextCand, nextCandESet, need - 2, newLastPos,
                     colorSeen, featSeen, stamp);
        }

        // 失活传播
        for (InstanceID w : activeNbrs) {
            {
                std::uint64_t k1 = edgeKey(u, w);
                auto it1 = _edgeId.find(k1);
                if (it1 != _edgeId.end()) {
                    markInactiveWithLog(_edgeAlive, it1->second, v, frameChangeLog);
                }
            }
            {
                std::uint64_t k2 = edgeKey(v, w);
                auto it2 = _edgeId.find(k2);
                if (it2 != _edgeId.end()) {
                    markInactiveWithLog(_edgeAlive, it2->second, u, frameChangeLog);
                }
            }
        }

        usedFeat[fu] = 0;
        usedFeat[fv] = 0;
        S.pop_back();
        S.pop_back();
    }

    // 恢复本次递归中修改的 alive 状态
    for (auto it = frameChangeLog.rbegin(); it != frameChangeLog.rend(); ++it) {
        _edgeAlive[it->first][it->second] = 1;
    }
}

std::uint64_t Model2::edgeKey(InstanceID a, InstanceID b) {
    if (a > b) std::swap(a, b);
    return (static_cast<std::uint64_t>(a) << 32) |
           static_cast<std::uint32_t>(b);
}

void Model2::execute() {
    buildColorMap();
    buildEdgeStructures();
    EdgeDecomposition();
    buildIncidentEdges();
    initEdgeAliveMasks();

    size = 0;
    std::vector<int> colorSeen(_features.size() + 1, 0);
    std::vector<int> featSeen(FeatureNumber.size() + 1, 0);
    int stamp = 0;

    std::vector<int> edgeSeen(_edges.size(), 0);
    int edgeSeenStamp = 0;

    for (int eid : _edgeOrder) {
        InstanceID u = _edges[eid].u;
        InstanceID v = _edges[eid].v;

        const int lastPos = _edgeRank[eid];

        // 永久删除：收集活跃公共邻居并标记相关边失效
        std::vector<InstanceID> seedActiveNeighbors;
        collectActiveNeighbors(eid, _edgeAlive, seedActiveNeighbors);
        for (InstanceID w : seedActiveNeighbors) {
            {
                std::uint64_t k1 = edgeKey(u, w);
                auto it1 = _edgeId.find(k1);
                if (it1 != _edgeId.end()) {
                    markInactive(_edgeAlive, it1->second, v);
                }
            }
            {
                std::uint64_t k2 = edgeKey(v, w);
                auto it2 = _edgeId.find(k2);
                if (it2 != _edgeId.end()) {
                    markInactive(_edgeAlive, it2->second, u);
                }
            }
        }

        // 候选顶点集
        std::vector<InstanceID> cand;
        collectActiveNeighbors(eid, _edgeAlive, cand);
        if ((int)cand.size() < _k - 2) continue;

        ++_vertexMarkStamp;
        const int vmark = _vertexMarkStamp;
        for (InstanceID x : cand) _vertexMark[x] = vmark;

        ++edgeSeenStamp;
        std::vector<int> canESet;
        canESet.reserve(cand.size() * 4);

        for (InstanceID x : cand) {
            for (int e2 : _incidentEdges[x]) {
                if (_edgeRank[e2] <= lastPos) continue;
                if (edgeSeen[e2] == edgeSeenStamp) continue;

                const InstanceID a = _edges[e2].u;
                const InstanceID b = _edges[e2].v;
                if (_vertexMark[a] != vmark || _vertexMark[b] != vmark) continue;
                if (_featureIdOfInstance[a] == _featureIdOfInstance[b]) continue;

                edgeSeen[e2] = edgeSeenStamp;
                canESet.push_back(e2);
            }
        }

        std::sort(canESet.begin(), canESet.end(),
                  [&](int e1, int e2) {
                      return _edgeRank[e1] < _edgeRank[e2];
                  });

        std::vector<InstanceID> S;
        S.reserve(_k);
        S.push_back(u);
        S.push_back(v);

        std::vector<uint8_t> usedFeat(FeatureNumber.size(), 0);
        usedFeat[_featureIdOfInstance[u]] = 1;
        usedFeat[_featureIdOfInstance[v]] = 1;

        // 直接使用 _edgeAlive，递归中会临时修改并通过日志恢复
        dfsEBBkC(S, usedFeat, cand, canESet, _k - 2, lastPos,
                 colorSeen, featSeen, stamp);
    }
}