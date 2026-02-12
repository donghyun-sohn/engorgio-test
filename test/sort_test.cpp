#include "ENGORGIO/comp.h"
#include "ENGORGIO/sort.h"
#include "ENGORGIO/quant.h"
#include "ENGORGIO/utils.h"
#include <stdexcept>
using namespace openfhe;
using namespace lbcrypto;

int main(int argc, char *argv[])
{
    int num_slots = 32768;  // default
    if (argc >= 2)
    {
        try { num_slots = std::stoi(argv[1]); }
        catch (...) { num_slots = 32768; }
    }
    std::cout << "Sorting " << num_slots << " records" << std::endl;

    std::chrono::system_clock::time_point start, end;
    start = std::chrono::system_clock::now();
    std::ifstream Sortfile("../../engorgio_sort.csv", std::ios::app);
    std::vector<std::vector<std::string>> csvData;
    std::string line;
    while (std::getline(Sortfile, line))
    {
        std::vector<std::string> row;
        std::stringstream ss(line);
        std::string cell;
        while (std::getline(ss, cell, ','))
        {
            row.push_back(cell);
        }
        csvData.push_back(row);
    }
    Sortfile.close();
    std::vector<double> sort_time;
    sort_time = bitonic_sort_modular(8, 2, num_slots);

    for (size_t i = 1; i < csvData.size(); i++)
    {
        csvData[i][5] = std::to_string(sort_time[i + 1]);
    }
    std::ofstream outFile_sort("../../engorgio_sort.csv");
    for (const auto &row : csvData)
    {
        for (size_t i = 0; i < row.size(); i++)
        {
            outFile_sort << row[i];
            if (i < row.size() - 1)
            {
                outFile_sort << ",";
            }
        }
        outFile_sort << std::endl;
    }
    outFile_sort.close();
    end = std::chrono::system_clock::now();
    double total_time = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << "Sort test compute time:" << total_time << std::endl;
    return 0;
}