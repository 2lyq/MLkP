#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

#include "CSVReader/CSVReader.h"
#include "EventClock/EventClock.hpp"
#include "Types.h"
#include "FindNeibs.h"
#include "Model2.h"
#include "FindRegion.h"

using namespace std;
namespace fs = filesystem;

int main(int argc, char **argv) {
    if (argc != 6) {
        cout << "Parameter format :\n"
             << argv[0]
             << " <dataset path> <minimum prevalent index> <candidate ratio> <distance> <k>\n";
        return 0;
    }

    string dataset(argv[1]);
    double minPre   = stod(argv[2]);
    double Ratio    = stod(argv[3]);
    double distance = stod(argv[4]);
    int    k        = stoi(argv[5]);

    InstanceMapper mapper;
    std::vector<InstanceID> instances;
    std::vector<LocationType> locations;
    std::vector<FeatureType> features;

    CSVReader csvReader(dataset, true);

    double min_x = std::numeric_limits<double>::max();
    double min_y = std::numeric_limits<double>::max();
    double max_x = -std::numeric_limits<double>::max();
    double max_y = -std::numeric_limits<double>::max();

    while (csvReader.hasNext()) {
        auto line = csvReader.getNextRecord();

        FeatureType feature = line[0];
        std::string idStr = line[1];
        double x = stod(line[2]);
        double y = stod(line[3]);
        LocationType location{x, y};

        InstanceID id = mapper.getID(feature, idStr);
        if (id >= static_cast<InstanceID>(locations.size())) {
            locations.resize(id + 1);
            features.resize(id + 1);
        }
        locations[id] = location;
        features[id] = feature;

        instances.push_back(id);

        if (x < min_x) min_x = x;
        if (y < min_y) min_y = y;
        if (x > max_x) max_x = x;
        if (y > max_y) max_y = y;
    }

    EventClock<TimeTicks::Microseconds> eventClock;
    eventClock.startClock("KRColocatingMining");

    FindNeibs findNeibs(instances, locations, features,
                        distance, min_x, min_y, max_x, max_y);
    NeibList &InsNeibs = findNeibs.InsNeibs;

    FindRegion findRegion(instances, InsNeibs, features, locations,
                          k, minPre, Ratio);

    auto RegionalPs = findRegion.execute();
    auto Colocation = findRegion.colocations;
    eventClock.stopClock("KRColocatingMining");
    const auto miningDuration = eventClock.getEventDuration("KRColocatingMining").count();

    double averageQS = 0.0;
    for (const auto &kv: RegionalPs) {
        averageQS += kv.second.second;
    }
    if (!RegionalPs.empty()) {
        averageQS /= static_cast<double>(RegionalPs.size());
    }

    {
        string filename = "PatternSummary" + to_string(distance) + ".txt";
        ofstream PatternSummary(filename);
        if (PatternSummary.is_open()) {
            PatternSummary << "PatternCount," << RegionalPs.size() << '\n';
            PatternSummary << "GlobalColocationCount," << Colocation.size() << '\n';
            PatternSummary << "AverageQS," << std::fixed << std::setprecision(10) << averageQS << '\n';
            PatternSummary << "RuntimeMicroseconds," << miningDuration << '\n';
        }
    }

    {
        ofstream Colocationtxt("Colocation.txt");
        if (Colocationtxt.is_open()) {
            for (const auto &pattern: Colocation) {
                for (const auto &f: pattern.feats) {
                    Colocationtxt << f;
                }
                Colocationtxt << '\n';
            }
        }
    }

    {
        string filename = "PatternsScore" + to_string(k) + ".txt";
        ofstream PatternsScore(filename);
        if (PatternsScore.is_open()) {
            for (const auto &kv: RegionalPs) {
                const ColocationType &pattern = kv.first;
                const auto &scores = kv.second.second;

                for (const auto &f: pattern.feats) {
                    PatternsScore << f;
                }
                PatternsScore << ',' << scores;
                PatternsScore << '\n';
            }
        }
    }

    {
        string filename = "MSRCount" + to_string(k) + ".txt";
        ofstream MSRCount(filename);
        if (MSRCount.is_open()) {
            std::size_t totalMSRCount = 0;
            for (const auto &kv: RegionalPs) {
                totalMSRCount += kv.second.first.size();
            }

            MSRCount << "TotalLCPPatterns," << RegionalPs.size() << '\n';
            MSRCount << "TotalMSRs," << totalMSRCount << '\n';
            MSRCount << "Pattern,MSRCount\n";

            for (const auto &kv: RegionalPs) {
                const ColocationType &pattern = kv.first;
                std::string patternName;
                for (const auto &f: pattern.feats) {
                    if (!patternName.empty()) patternName += "_";
                    patternName += f;
                }
                MSRCount << patternName << ',' << kv.second.first.size() << '\n';
            }
        }
    }

    fs::path outDir = fs::path("Regional Colocation_") / std::to_string(k);
    fs::create_directories(outDir);

    for (const auto &r: RegionalPs) {
        const ColocationType &pattern = r.first;
        const auto &regionList = r.second.first;

        std::string fileName;
        for (const auto &f: pattern.feats) {
            if (!fileName.empty()) fileName += "_";
            fileName += f;
        }

        for (std::size_t regionIndex = 0; regionIndex < regionList.size(); ++regionIndex) {
            const Region &regionIns = regionList[regionIndex];
            std::string outputName = fileName;
            if (regionList.size() > 1) {
                outputName += "_MSR" + std::to_string(regionIndex + 1);
            }

            fs::path filePath = outDir / (outputName + ".csv");
            std::ofstream Regions(filePath);
            if (!Regions.is_open()) {
                std::cerr << "Unable to create file: " << filePath.string() << std::endl;
                continue;
            }

            Regions << "Feature,Instance,LocationX,LocationY\n";
            for (auto id: regionIns) {
                const auto &feature = features[id];
                const auto &label = mapper.getLabel(id);
                const auto &loc = locations[id];

                Regions << feature << ","
                        << label << ","
                        << std::fixed << std::setprecision(6) << loc.first << ","
                        << std::fixed << std::setprecision(6) << loc.second << "\n";
            }
        }
    }

    return 0;
}
