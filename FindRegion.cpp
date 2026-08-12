#include "FindRegion.h"
#include <numeric>
#include <queue>
#include <cmath>
#include <algorithm>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <functional>
#include <climits>

// -----------------------------------------------------
// 构造函数
// -----------------------------------------------------
FindRegion::FindRegion(const std::vector<InstanceID>& instances,
                       const NeibList& adj,
                       const std::vector<FeatureType>& features,
                       const std::vector<LocationType>& locations,
                       int k,
                       double threshold,
                       double ratio)
        : _k(k),
          _Threshold(threshold),
          _Ratio(ratio),
          _instances(instances),
          _adj(adj),
          _features(&features),
          _locations(locations),
          _model(instances, adj, features, k)
{
    for (InstanceID id : _instances) {
        if (id < 0 || static_cast<std::size_t>(id) >= features.size()) continue;
        const FeatureType& f = features[id];
        ++FeatureNumber[f];
    }
    setupPatternCollector();
    _model.execute();
}


// -----------------------------------------------------
void FindRegion::setupPatternCollector() {
    _model.onClique = [this](const std::vector<InstanceID>& S) {
        // 防止 _features 为空
        if (!this->_features) return;

        const auto& featsById = *this->_features;
        const std::size_t n = featsById.size();

        std::vector<FeatureType> feats;
        feats.reserve(S.size());

        for (auto id : S) {
            if (id < 0) continue;
            const std::size_t uid = static_cast<std::size_t>(id);
            if (uid >= n) continue;

            feats.push_back(featsById[uid]);
        }

        // 如果过滤后不足以构成模式，可直接退出（可选）
        if (feats.empty()) return;

        ColocationType pattern(feats);
        auto& W = PatternW[pattern];

        for (auto id : S) {
            if (id < 0) continue;
            const std::size_t uid = static_cast<std::size_t>(id);
            if (uid >= n) continue;

            W[id] += 1.0;
        }
    };
}

// -----------------------------------------------------
FindRegion::WeightMapType
FindRegion::SetWeight(const ColocationType& pattern,
                      const WeightMapType& rawW) const {
    WeightMapType W;
    for (auto id : _instances) {
        if (std::find(pattern.feats.begin(), pattern.feats.end(),
                      (*_features)[id]) != pattern.feats.end()) {
            W[id] = (rawW.count(id) ? rawW.at(id) : 0.0);

        }
    }
    return W;
}

// -----------------------------------------------------
bool FindRegion::isGColocation(const WeightMapType& Weight) const {
    std::unordered_map<FeatureType,std::unordered_set<InstanceID>> Fs;
    const auto& feats = *_features;
    const std::size_t n = feats.size();
    for (const auto& kv : Weight) {
        if (kv.second <= 0) continue;

        InstanceID id = kv.first;
        if (id < 0) continue;

        const std::size_t uid = static_cast<std::size_t>(id);
        if (uid >= n) continue;

        Fs[feats[uid]].insert(id);
    }

    for (auto& kv : Fs) {
        double ratio = kv.second.size() * 1.0 / FeatureNumber.at(kv.first);
        if (ratio < _Threshold) return false;
    }
    return true;
}

// -----------------------------------------------------
bool FindRegion::isCandidate(const WeightMapType& Weight,
                             const ColocationType& pattern) const {
    std::unordered_map<FeatureType,int> cnt;
    const auto& feats = *_features;
    const std::size_t n = feats.size();

    for (const auto& kv : Weight) {
        if (kv.second <= 0) continue;

        InstanceID id = kv.first;
        if (id < 0) continue;

        const std::size_t uid = static_cast<std::size_t>(id);
        if (uid >= n) continue;

        ++cnt[feats[uid]];
    }
    for (auto& f : pattern.feats) {
        double ratio = cnt[f] * 1.0 / FeatureNumber.at(f);
        if (ratio < _Ratio) return false;
    }
    return true;
}


// 统计参与实例最少的特征作为种子特征
FeatureType FindRegion::selectSeedFeature(const InstanceWeightMap<double>& WeightMap,
                                          const ColocationType& Candidate)
{
    // 1) 将模式特征放到一个集合里，快速判断 “是否属于该模式”
    std::unordered_set<FeatureType> candSet(Candidate.feats.begin(),
                                            Candidate.feats.end());

    // 2) 统计每个特征的参与实例数量
    std::unordered_map<FeatureType,int> FNumber;

    for (const auto& kv : WeightMap) {
        InstanceID ins = kv.first;
        double     w   = kv.second;
        if (w <= 0.0) continue;  // 只统计参与实例

        if (ins < 0 || static_cast<std::size_t>(ins) >= _features->size()) continue;

        const FeatureType& f = (*_features)[ins];

        // 只统计模式 C 中的特征
        if (candSet.find(f) != candSet.end()) {
            ++FNumber[f];
        }
    }

    // 3) 在模式特征中选择参与实例数量最少的那个
    FeatureType seedFeature = Candidate.feats[0];
    int minSize = INT_MAX;

    for (const auto& f : Candidate.feats) {
        int cnt = 0;
        auto it = FNumber.find(f);
        if (it != FNumber.end()) {
            cnt = it->second;
        }
        if (cnt < minSize) {
            minSize    = cnt;
            seedFeature = f;
        }
    }

    return seedFeature;
}
RegionStats FindRegion::buildRegionStatsFromScratch(
        const Region &reg,
        const std::vector<double> &weight,
        const ColocationType &Candidate,
        const std::unordered_map<FeatureType,int> &feat2idx,
        const std::vector<int> &globalCnt,
        double N_C)
{
    RegionStats st;
    const int k = static_cast<int>(Candidate.feats.size());
    st.cntR.assign(k, 0);
    st.cntP.assign(k, 0);

    st.N_R = static_cast<int>(reg.size());

    // 1) sumW, cntR, cntP
    for (InstanceID ins : reg) {
        if (ins < 0 || static_cast<std::size_t>(ins) >= _features->size()) continue;

        const FeatureType &f = (*_features)[ins];
        auto itIdx = feat2idx.find(f);
        if (itIdx == feat2idx.end()) continue;

        int idx = itIdx->second;
        ++st.cntR[idx];

        if (ins >= 0 && ins < (int)weight.size() && weight[ins] > 0.0) {
            ++st.cntP[idx];
            st.sumW += weight[ins];
        }
    }

    // 2) RI, N_CR, RPI + ?
    st.RI   = st.sumW / std::max(1, _k);
    st.N_CR = 0.0;
    st.RPI  = 1.0;

    const int ksize = static_cast<int>(Candidate.feats.size());
    for (int i = 0; i < ksize; ++i) {
        int nr = st.cntR[i];
        int gc = globalCnt[i];

        st.N_CR += nr;

        double RPIf = 0.0;
        if (nr > 0) {
            RPIf = static_cast<double>(st.cntP[i]) / static_cast<double>(nr);
        }

        if (RPIf < st.RPI) st.RPI = RPIf;
    }

    if (st.RI <= 0.0 || st.N_R <= 0 || N_C <= 0.0) {
        st.score = 0.0;
        return st;
    }

    long double logv  = st.RI/allRI;

    long double gamma = 1.0L / (long double)_k;
    double term1_2= (std::pow(st.RI,gamma))/st.N_CR;
    term1_2 = std::pow(term1_2,gamma);
    st.score = (logv+term1_2) / 2.0;
    return st;
}
RegionStats FindRegion::mergeStatsIncremental(
        const RegionStats &a,
        const RegionStats &b,
        const std::vector<InstanceID> &bridge,
        const std::vector<double> &weight,
        const ColocationType &Candidate,
        const std::unordered_map<FeatureType,int> &feat2idx,
        const std::vector<int> &globalCnt,
        double N_C)
{
    RegionStats st;
    const int k = static_cast<int>(Candidate.feats.size());
    st.cntR.assign(k, 0);
    st.cntP.assign(k, 0);

    // 合并两个区域块的统计量
    st.sumW = a.sumW + b.sumW;
    st.N_R  = a.N_R  + b.N_R;
    st.N_CR = a.N_CR + b.N_CR;

    for (int i = 0; i < k; ++i) {
        st.cntR[i] = a.cntR[i] + b.cntR[i];
        st.cntP[i] = a.cntP[i] + b.cntP[i];
    }

    // 补充两个区域块之间的桥接实例
    for (InstanceID x : bridge) {
        if (x < 0 || x >= (int)weight.size()) continue;
        if (static_cast<std::size_t>(x) >= _features->size()) continue;

        const FeatureType &f = (*_features)[x];
        auto itIdx = feat2idx.find(f);
        if (itIdx == feat2idx.end()) continue;

        int idx = itIdx->second;
        ++st.cntR[idx];
        ++st.N_R;
        st.N_CR += 1;

        if (weight[x] > 0.0) {
            ++st.cntP[idx];
            st.sumW += weight[x];
        }
    }

    // 重新计算合并区域的 RI 和 RPI
    st.RI  = st.sumW / std::max(1, _k);
    st.RPI = 1.0;

    for (int i = 0; i < k; ++i) {
        int nr = st.cntR[i];
        int gc = globalCnt[i];

        double RPIf = 0.0;
        if (nr > 0) {
            RPIf = static_cast<double>(st.cntP[i]) / static_cast<double>(nr);
        }

        if (RPIf < st.RPI) st.RPI = RPIf;
    }

    if (st.RI <= 0.0 || st.N_R <= 0 || N_C <= 0.0) {
        st.score = 0.0;
        return st;
    }

    long double logv  = st.RI/allRI;
    long double gamma = 1.0L / (long double)_k;

    double term1_2= (std::pow(st.RI,gamma))/st.N_CR;
    term1_2 = std::pow(term1_2,gamma);
    st.score = (logv+term1_2 ) / 2.0;
    return st;
}




std::pair<FindRegion::RegionList,double>
FindRegion::identify_best_region(const InstanceWeightMap<double>& W,
                                 const ColocationType&            Candidate,
                                 const FeatureType&               seedFeature)
{
    int n = static_cast<int>(_adj.size());
    if (n <= 0) return {{}, 0.0};

    // 0) weight  unordered_map
    std::vector<double> weight(n, 0.0);
    for (const auto &kv : W) {
        InstanceID ins = kv.first;
        double     ww  = kv.second;
        if (ins < 0 || ins >= n) continue;
        if (ww > 0.0) weight[ins] = ww;
    }

    std::vector<char> isParticipating(n, false);
    for (int i = 0; i < n; ++i) {
        if (weight[i] > 0.0) isParticipating[i] = true;
    }

    // 1) 选择种子特征实例构建block
    std::vector<InstanceID> seeds;
    if (_features) {
        seeds.reserve(64);
        for (const auto &kv : W) {
            InstanceID ins = kv.first;
            if (ins < 0 || ins >= n) continue;
            if (weight[ins] <= 0.0) continue;
            if (static_cast<std::size_t>(ins) >= _features->size()) continue;
            if ((*_features)[ins] == seedFeature) {
                seeds.push_back(ins);
            }
        }
    }

    if (seeds.empty()) {
        return {{}, 0.0};
    }

    // 2) BFS to get block
    std::vector<char> visited(n, false);
    std::vector<Region> regions;
    regions.reserve(seeds.size());

    std::queue<InstanceID> q;
    for (InstanceID s : seeds) {
        if (s < 0 || s >= n) continue;
        if (!isParticipating[s]) continue;
        if (visited[s]) continue;

        Region reg;
        reg.reserve(128);

        visited[s] = true;
        q.push(s);
        reg.push_back(s);

        while (!q.empty()) {
            InstanceID u = q.front(); q.pop();
            if (u < 0 || u >= n) continue;

            const auto &nbrs = _adj[u];
            for (InstanceID nb : nbrs) {
                if (nb < 0 || nb >= n) continue;
                if (!isParticipating[nb]) continue;
                if (visited[nb]) continue;
                visited[nb] = true;
                q.push(nb);
                reg.push_back(nb);
            }
        }

        if (!reg.empty()) {
            regions.push_back(std::move(reg));
        }
    }

    int m = static_cast<int>(regions.size());
    if (m == 0) return {{}, 0.0};

    // 构建模式特征索引和全局实例数量
    const int k = static_cast<int>(Candidate.feats.size());
    if (k <= 0) return {{}, 0.0};

    std::unordered_map<FeatureType,int> feat2idx;
    feat2idx.reserve(k * 2);
    for (int i = 0; i < k; ++i) {
        feat2idx[Candidate.feats[i]] = i;
    }

    std::vector<int> globalCnt(k, 0);
    double N_C = 0.0;
    for (int i = 0; i < k; ++i) {
        const auto &f = Candidate.feats[i];
        auto itG = FeatureNumber.find(f);
        if (itG != FeatureNumber.end()) {
            globalCnt[i] = itG->second;
            N_C += itG->second;
        }
    }
    if (N_C <= 0.0) return {{}, 0.0};

    // ============================================================
    // 计算每个block的值
    // ============================================================

    std::vector<RegionStats> stats(m);
    std::vector<double>      sumW(m, 0.0);
    std::vector<double>      score(m, 0.0);

    // 2.1 getRegionQualityScore
    for (int i = 0; i < m; ++i) {
        stats[i] = buildRegionStatsFromScratch(
                regions[i],
                weight,
                Candidate,
                feat2idx,
                globalCnt,
                N_C);

        sumW[i] = stats[i].sumW;
        score[i] = stats[i].score;
    }

    // 2.2 instance -> regionID
    std::unordered_map<InstanceID,int> insToRegion;
    insToRegion.reserve(1024);
    for (int rid = 0; rid < m; ++rid) {
        for (InstanceID ins : regions[rid]) {
            insToRegion[ins] = rid;
        }
    }

    // 2.3 3跳获取邻居block
    std::vector<std::unordered_set<int>> regionNbrs(m);
    using BridgeMap = std::unordered_map<int, std::vector<InstanceID>>;
    std::vector<BridgeMap> regionBridge(m);

    for (int rid = 0; rid < m; ++rid) {
        for (InstanceID oa : regions[rid]) {
            if (oa < 0 || oa >= n) continue;

            const auto &nbr_oa = _adj[oa];
            for (InstanceID oi : nbr_oa) {
                if (oi < 0 || oi >= n) continue;

                const auto &nbr_oi = _adj[oi];
                for (InstanceID ij : nbr_oi) {
                    if (ij < 0 || ij >= n) continue;

                    const auto &nbr_ij = _adj[ij];
                    for (InstanceID ob : nbr_ij) {
                        if (ob < 0 || ob >= n) continue;

                        auto it = insToRegion.find(ob);
                        if (it == insToRegion.end()) continue;
                        int rid2 = it->second;
                        if (rid2 == rid) continue;

                        regionNbrs[rid].insert(rid2);
                        regionNbrs[rid2].insert(rid);

                        auto &bm1 = regionBridge[rid];
                        if (bm1.find(rid2) == bm1.end()) {
                            bm1[rid2] = {oi, ij};
                        }
                        auto &bm2 = regionBridge[rid2];
                        if (bm2.find(rid) == bm2.end()) {
                            bm2[rid] = {oi, ij};
                        }
                    }
                }
            }
        }
    }



// ============================================================
// 2.4 Branch-and-Bound (B&B) search on region blocks
//     - Each initial connected component (regions[i]) is treated as a block
//     - Objective: maximize two-term RCS = (RIC + SRI) / 2
//       where RIC = RI / allRI, and SRI is the compressed closeness term.
//     - Constraint: RPI >= _Threshold (RPI is NOT part of the objective)
// ============================================================

// Convert neighbor sets to adjacency lists (block graph)
    std::vector<std::vector<int>> blockAdj(m);
    blockAdj.reserve(m);
    for (int i = 0; i < m; ++i) {
        blockAdj[i].assign(regionNbrs[i].begin(), regionNbrs[i].end());
    }

// Block potential: sumW is the sum of participating-instance weights in the block
    std::vector<double> blockW(m, 0.0);
    for (int i = 0; i < m; ++i) {
        blockW[i] = stats[i].sumW;
    }

// Rank blocks by their current block-level region score (descending).
// This is a heuristic ordering to find good solutions early; it does not change correctness.
    std::vector<int> blockOrder(m);
    std::iota(blockOrder.begin(), blockOrder.end(), 0);
    std::stable_sort(blockOrder.begin(), blockOrder.end(),
                     [&](int a, int b) {
                         if (stats[a].score != stats[b].score) return stats[a].score > stats[b].score;
                         return a < b;
                     });

    std::vector<int> rank(m, 0);
    for (int pos = 0; pos < m; ++pos) rank[blockOrder[pos]] = pos;

// Symmetry breaking on score-rank:
// when root is r, only blocks with rank >= rank[r] are allowed.
// Hence each connected block-subgraph is explored exactly once (at its highest-score block).
    auto allowed = [&](int b, int root) -> bool {
        return rank[b] >= rank[root];
    };

// Helper: compute reachable "remaining potential" (tight upper bound under connectivity + excluded + rank constraint).
// We compute the connected component of root in the induced block graph over allowed & non-excluded blocks.
// Any future connected expansion must lie within this component.
    auto computeReachableRemaining = [&](int root,
                                         const std::vector<char>& chosen,
                                         const std::vector<char>& excluded) -> std::pair<double, std::vector<int>>
    {
        std::vector<char> seen(m, false);
        std::queue<int> qq;
        seen[root] = true;
        qq.push(root);

        while (!qq.empty()) {
            int u = qq.front(); qq.pop();
            for (int nb : blockAdj[u]) {
                if (!allowed(nb, root)) continue;
                if (excluded[nb]) continue;
                if (seen[nb]) continue;
                seen[nb] = true;
                qq.push(nb);
            }
        }

        double remW = 0.0;
        std::vector<int> remP(_k, 0);
        for (int b = 0; b < m; ++b) {
            if (!seen[b]) continue;
            if (!allowed(b, root)) continue;
            if (chosen[b] || excluded[b]) continue;

            remW += blockW[b];
            for (int j = 0; j < _k; ++j) remP[j] += stats[b].cntP[j];
        }
        return {remW, remP};
    };

// Upper bound for RPI feasibility (constraint).
// Freeze denominators cntR, and optimistically add all remaining participating counts to numerators.
    auto computeRPIUpperBound = [&](const RegionStats& cur,
                                    const std::vector<int>& remP) -> double
    {
        long double rpiUB = 1.0L;
        for (int i = 0; i < _k; ++i) {
            const int denom = cur.cntR[i];
            long double ub = 0.0L;

            if (denom > 0) {
                long double numer = (long double)cur.cntP[i] + (long double)remP[i];
                ub = numer / (long double)denom;
                if (ub > 1.0L) ub = 1.0L;
                if (ub < 0.0L) ub = 0.0L;
            } else {
                // no instance of this feature in current region;
                // if any remaining participating exists, optimistic UB=1, otherwise 0
                ub = (remP[i] > 0) ? 1.0L : 0.0L;
            }

            if (ub < rpiUB) rpiUB = ub;
        }
        if (rpiUB < 0.0L) rpiUB = 0.0L;
        if (rpiUB > 1.0L) rpiUB = 1.0L;
        return (double)rpiUB;
    };

// Upper bound on score (objective): optimistic RI -> optimistic RIC and SRI.
// - For RI upper bound, assume all remaining block weights can be added.
// - For SRI upper bound, maximize RI and freeze N_CR to current value (expansion only increases N_CR, so this is optimistic).
    auto computeScoreUpperBound = [&](const RegionStats& cur, double remW) -> double
    {
        if (allRI <= 0.0L) return 0.0;

        const long double gamma = 1.0L / (long double)_k;

        // optimistic RI upper bound
        long double RIub = ((long double)cur.sumW + (long double)remW) / (long double)std::max(1, _k);
        if (RIub < 0.0L) RIub = 0.0L;

        // RIC upper bound: RI / allRI
        long double RICub = RIub / (long double)allRI;
        if (RICub < 0.0L) RICub = 0.0L;
        if (RICub > 1.0L) RICub = 1.0L;

        // SRI upper bound: ((RI^{1/k}) / N_CR)^{1/k}, with N_CR frozen to current (optimistic)
        long double denom = (long double)cur.N_CR;
        if (denom < 1.0L) denom = 1.0L;

        long double tmp = 0.0L;
        if (RIub > 0.0L) {
            tmp = std::pow(RIub, gamma) / denom;
            if (tmp < 0.0L) tmp = 0.0L;
        }

        long double SRIub = std::pow(tmp, gamma);
        if (SRIub < 0.0L) SRIub = 0.0L;
        if (SRIub > 1.0L) SRIub = 1.0L;

        long double ubScore = (RICub + SRIub) / 2.0L;
        if (ubScore < 0.0L) ubScore = 0.0L;
        if (ubScore > 1.0L) ubScore = 1.0L;
        return (double)ubScore;
    };

// Global best (across all roots)
    constexpr double ScoreEps = 1e-12;
    double bestScore = 0.0;
    struct BestSolution {
        std::vector<int> blocks;
        std::vector<int> parent;
    };
    std::vector<BestSolution> bestSolutions;

// Depth-first Branch-and-Bound on the block graph (include/exclude a frontier block)
    std::function<void(
            int,                      // root block id
            std::vector<char>&,        // chosen blocks
            std::vector<char>&,        // excluded blocks
            std::vector<char>&,        // inFrontier
            std::vector<int>&,         // frontier list (may contain stale entries; guarded by inFrontier)
            RegionStats&,              // current aggregated stats
            std::vector<int>&          // parent block for reconstruction
    )> dfs;

    dfs = [&](int root,
              std::vector<char>& chosen,
              std::vector<char>& excluded,
              std::vector<char>& inFrontier,
              std::vector<int>& frontier,
              RegionStats& cur,
              std::vector<int>& parent)
    {
        // Update best feasible solution (constraint satisfied)
        if (cur.RPI >= _Threshold && cur.score > 0.0) {
            if (cur.score > bestScore + ScoreEps) {
                bestScore = cur.score;
                bestSolutions.clear();
            }
            if (std::fabs(cur.score - bestScore) <= ScoreEps) {
                BestSolution sol;
                sol.parent = parent;
                for (int i = 0; i < m; ++i) if (chosen[i]) sol.blocks.push_back(i);
                bestSolutions.push_back(std::move(sol));
            }
        }

        // Tight remaining potentials under current excluded + rank + connectivity
        auto rem = computeReachableRemaining(root, chosen, excluded);
        const double remW = rem.first;
        const std::vector<int>& remP = rem.second;

        // Feasibility pruning (RPI cannot reach the threshold)
        const double rpiUB = computeRPIUpperBound(cur, remP);
        if (rpiUB < _Threshold - 1e-12) return;

        // Score upper bound pruning
        const double ubScore = computeScoreUpperBound(cur, remW);
        if (ubScore < bestScore - ScoreEps) return;

        // pick a frontier block (heuristic: max block-level score)
        int pick = -1;
        double pickScore = -1.0;
        for (int b : frontier) {
            if (!inFrontier[b]) continue;
            if (chosen[b] || excluded[b]) continue;
            if (!allowed(b, root)) continue;

            if (stats[b].score > pickScore) {
                pickScore = stats[b].score;
                pick = b;
            }
        }
        if (pick == -1) return;

        // remove pick from frontier for branching (lazy deletion via inFrontier)
        inFrontier[pick] = false;

        // ----------------------
        // Branch 1: INCLUDE pick
        // ----------------------
        {
            std::vector<char> chosen2   = chosen;
            std::vector<char> excluded2 = excluded;
            std::vector<char> inF2      = inFrontier;
            std::vector<int>  frontier2 = frontier;
            std::vector<int>  parent2   = parent;
            RegionStats        cur2     = cur;

            chosen2[pick] = true;

            // find a connector block already chosen and adjacent to pick
            int conn = -1;
            for (int nb : blockAdj[pick]) {
                if (!allowed(nb, root)) continue;
                if (chosen[nb]) { conn = nb; break; }
            }

            if (conn != -1) {
                parent2[pick] = conn;

                // bridge instances that connect conn -> pick (length-2 bridge stored in regionBridge)
                std::vector<InstanceID> bridge;
                auto itBr = regionBridge[conn].find(pick);
                if (itBr != regionBridge[conn].end()) bridge = itBr->second;

                cur2 = mergeStatsIncremental(cur, stats[pick], bridge,
                                             weight, Candidate, feat2idx, globalCnt, N_C);

                // update frontier: add neighbors of pick
                for (int nb : blockAdj[pick]) {
                    if (!allowed(nb, root)) continue;
                    if (chosen2[nb] || excluded2[nb]) continue;
                    if (!inF2[nb]) {
                        inF2[nb] = true;
                        frontier2.push_back(nb);
                    }
                }

                dfs(root, chosen2, excluded2, inF2, frontier2, cur2, parent2);
            }
        }

        // ----------------------
        // Branch 2: EXCLUDE pick
        // ----------------------
        {
            excluded[pick] = true;
            dfs(root, chosen, excluded, inFrontier, frontier, cur, parent);
        }
    };

// Run B&B from each block as root, in descending score order (blockOrder).
    for (int t = 0; t < m; ++t) {
        int root = blockOrder[t];

        std::vector<char> chosen(m, false), excluded(m, false), inFrontier(m, false);
        std::vector<int>  parent(m, -1);

        chosen[root] = true;

        // Initial frontier: neighbors of root allowed by rank constraint
        std::vector<int> frontier;
        frontier.reserve(blockAdj[root].size());
        for (int nb : blockAdj[root]) {
            if (!allowed(nb, root)) continue;
            if (!inFrontier[nb] && !chosen[nb] && !excluded[nb]) {
                inFrontier[nb] = true;
                frontier.push_back(nb);
            }
        }

        RegionStats cur = stats[root];
        dfs(root, chosen, excluded, inFrontier, frontier, cur, parent);
    }
    if (bestSolutions.empty()) {
        return {{}, 0.0};
    }

// Reconstruct all best region instances = union of selected blocks + bridge instances on selected edges
    RegionList finalRegions;
    finalRegions.reserve(bestSolutions.size());

    for (const auto& sol : bestSolutions) {
        Region finalReg;
        finalReg.reserve(256);
        std::vector<char> inIns(n, false);

        auto addIns = [&](InstanceID id) {
            if (id < 0 || id >= n) return;
            if (!inIns[id]) {
                inIns[id] = true;
                finalReg.push_back(id);
            }
        };

        for (int bid : sol.blocks) {
            for (InstanceID id : regions[bid]) addIns(id);
        }
        for (int bid : sol.blocks) {
            int par = sol.parent[bid];
            if (par >= 0) {
                auto itBr = regionBridge[par].find(bid);
                if (itBr != regionBridge[par].end()) {
                    for (InstanceID id : itBr->second) addIns(id);
                }
            }
        }

        finalRegions.push_back(std::move(finalReg));
    }

    return {finalRegions, bestScore};

}



std::map<ColocationType,std::pair<FindRegion::RegionList,double>>
FindRegion::execute()
{
    std::map<ColocationType,std::pair<RegionList,double>> result;

    for (auto& kv : PatternW) {
        const ColocationType& pattern = kv.first;
        const WeightMapType& rawW = kv.second;
        allRI = 0;
        auto W = SetWeight(pattern, rawW);
        for(auto ins:W){
            allRI+=ins.second;
        }
        allRI=allRI*1.0/_k;
        if (isGColocation(W)) {
            colocations.insert(pattern);
            continue;
        }

        if (!isCandidate(W, pattern)) continue;

        FeatureType seedFeature = selectSeedFeature(W,pattern);
        auto best = identify_best_region(W,pattern,seedFeature);

        if (best.second > 0)
            result[pattern] = best;
    }

    return result;
}

