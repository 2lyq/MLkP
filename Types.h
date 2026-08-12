#pragma once
#include <string>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <cstdint>
#include <algorithm>

/* ---------------------------------------------------------
   1. 基础类型
   --------------------------------------------------------- */

using InstanceID = int;                  // 实例 ID（统一用 int）
using FeatureType = std::string;         // 特征类型：string
using LocationType = std::pair<double,double>;

using CellId = std::pair<int,int>;
using Region = std::vector<InstanceID>;

/* ---------------------------------------------------------
   2. ColocationType（模式类型）
   --------------------------------------------------------- */

struct ColocationType {
    std::vector<FeatureType> feats;

    ColocationType() = default;

    explicit ColocationType(const std::vector<FeatureType> &v) {
        feats = v;
        std::sort(feats.begin(), feats.end());
        feats.erase(std::unique(feats.begin(), feats.end()), feats.end());
    }

    bool operator==(const ColocationType &other) const {
        return feats == other.feats;
    }

    bool operator<(const ColocationType &other) const {
        return feats < other.feats;
    }
};

struct ColocationHasher {
    size_t operator()(ColocationType const &c) const noexcept {
        std::size_t h = 0;
        for (auto &s : c.feats) {
            h ^= std::hash<std::string>{}(s)
                 + 0x9e3779b97f4a7c15ULL
                 + (h << 6)
                 + (h >> 2);
        }
        return h;
    }
};

/* ---------------------------------------------------------
   3. InstanceMapper：string ↔ InstanceID 映射
   --------------------------------------------------------- */
class InstanceMapper {
public:
    InstanceID getID(const std::string &feature, const std::string &label) {
        std::string key = feature + "#" + label;
        auto it = name2id.find(key);
        if (it != name2id.end()) return it->second;

        InstanceID newId = (InstanceID)id2feature.size();
        name2id[key] = newId;

        id2feature.push_back(feature);
        id2label.push_back(label);
        return newId;
    }

    const std::string &getFeature(InstanceID id) const {
        return id2feature[id];
    }

    const std::string &getLabel(InstanceID id) const {
        return id2label[id];
    }

    int size() const { return (int)id2feature.size(); }

private:
    std::unordered_map<std::string, InstanceID> name2id;
    std::vector<std::string> id2feature;
    std::vector<std::string> id2label;
};

/* ---------------------------------------------------------
   4. 邻接表
   --------------------------------------------------------- */

using NeibList = std::vector<std::vector<InstanceID>>;

/* ---------------------------------------------------------
   5. Pattern → (InstanceID → 权重)
   --------------------------------------------------------- */

template<typename T>
using InstanceWeightMap = std::unordered_map<InstanceID, T>;

using PatternWType = std::unordered_map<
        ColocationType,
        InstanceWeightMap<double>,
        ColocationHasher>;

/* ---------------------------------------------------------
   6. Color 类型
   --------------------------------------------------------- */
using Color = int;
struct RegionStats {
    int    N_R   = 0;     // |R|
    double sumW  = 0.0;   // 区域总权重（只计 W>0）
    double RI    = 0.0;   // 行实例（近似 sumW / k）
    double N_CR  = 0.0;   // sum_f N_f(R)
    double RPI   = 0.0;   // min_f Np_f / N_f
    double score = 0.0;   // RCS(C,R)

    std::vector<int> cntR; // N_f(R) size = k
    std::vector<int> cntP; // Np_f(R) size = k
};

