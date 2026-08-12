#pragma once

#include <vector>
#include <map>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <utility>

#include "Types.h"
#include "Model2.h"

// 区域模式挖掘
class FindRegion {
public:
    using WeightMapType = InstanceWeightMap<double>;   // InstanceID -> weight
    using Region = std::vector<InstanceID>;
    using RegionList = std::vector<Region>;

    FindRegion(const std::vector<InstanceID>& instances,
               const NeibList& adj,
               const std::vector<FeatureType>& features,
               const std::vector<LocationType>& locations,
               int k,
               double threshold,
               double ratio);

    // 主流程：返回 (pattern -> {所有最大 QS 区域, Sr})
    std::map<ColocationType, std::pair<RegionList,double>> execute();

    const std::map<FeatureType,int>& getFeatureNumber() const { return FeatureNumber; }
    const std::map<ColocationType, WeightMapType>& getPatternW() const { return PatternW; }
    const std::map<ColocationType, std::vector<double>>& getPatternScore() const { return PatternScore; }
    const std::set<ColocationType>& getGlobalColocations() const { return colocations; }

private:
    // 构造 PatternW：模式实例权重
    void buildFeatureNumber();
    void setupPatternCollector();

    WeightMapType SetWeight(const ColocationType& pattern,
                            const WeightMapType& rawW) const;

    bool isGColocation(const WeightMapType& Weight) const;
    bool isCandidate(const WeightMapType& Weight,
                     const ColocationType& pattern) const;

    // 扩展候选区域
    std::vector<Region> grow_regions(const WeightMapType& W,
                                     const FeatureType& seedFeature) const;

    std::vector<double> Calculate(const std::pair<Region,double>& RegionW,
                                  const WeightMapType& W,
                                  const ColocationType& Candidate);

    std::pair<Region,double> GetOptimalRegion(const std::vector<Region>& regions,
                                              const WeightMapType& W,
                                              const ColocationType& Candidate);
    FeatureType selectSeedFeature(const InstanceWeightMap<double>& WeightMap,
                                  const ColocationType& Candidate);
    std::pair<RegionList,double> identify_best_region(const InstanceWeightMap<double>& W,
                                                       const ColocationType&            Candidate,
                                                       const FeatureType&               seedFeature);
    RegionStats buildRegionStatsFromScratch(
            const Region &reg,
            const std::vector<double> &weight,
            const ColocationType &Candidate,
            const std::unordered_map<FeatureType,int> &feat2idx,
            const std::vector<int> &globalCnt,
            double N_C);
    RegionStats mergeStatsIncremental(
            const RegionStats &a,
            const RegionStats &b,
            const std::vector<InstanceID> &bridge,
            const std::vector<double> &weight,
            const ColocationType &Candidate,
            const std::unordered_map<FeatureType,int> &feat2idx,
            const std::vector<int> &globalCnt,
            double N_C);

public:
    int    _k;
    double _Threshold;
    double _Ratio;
    int allRI;

    // 全局实例信息
    std::vector<InstanceID>   _instances;
    NeibList                  _adj;
    const std::vector<FeatureType>*  _features;
    std::vector<LocationType> _locations;

    // 特征全局数量
    std::map<FeatureType,int> FeatureNumber;

    // 模式 -> (InstanceID -> weight)
    std::map<ColocationType, WeightMapType> PatternW;

    // 模式评分：Sr、Sr_1、RPI、Rn/Mu2
    std::map<ColocationType, std::vector<double>> PatternScore;

    // 全局 G-colocations
    std::set<ColocationType> colocations;

    // 内部模型
    Model2 _model;
};

