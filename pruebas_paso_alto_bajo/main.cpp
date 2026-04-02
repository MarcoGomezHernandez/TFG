#include <fstream>
#include <iostream>
#include <string>
#include <algorithm>
#include <kfr/all.hpp>

int main(int argc, char *argv[])
{
    if (argc < 6)
        return std::cerr << "Uso: " << argv[0] << " <file.csv> <fs> <fc> <npoints> <col> [--skip-header]\n", 1;

    double fs, fc;
    int npoints, col;
    try
    {
        fs = std::stod(argv[2]);
        fc = std::stod(argv[3]);
        npoints = std::stoi(argv[4]);
        col = std::stoi(argv[5]);
    }
    catch (...)
    {
        return 1;
    }

    std::ifstream file(argv[1]);
    std::string line;
    if (argc > 6)
        std::getline(file, line);

    kfr::univector<double> data;
    while (std::getline(file, line) && data.size() < (size_t)npoints)
    {
        size_t pos = 0;
        for (int i = 0; i < col && pos != std::string::npos; ++i)
        {
            pos = line.find(',', pos);
            if (pos != std::string::npos)
                pos++;
        }
        if (pos != std::string::npos)
        {
            try
            {
                data.push_back(std::stod(line.substr(pos)));
            }
            catch (...)
            {
            }
        }
    }

    if (data.empty())
        return 1;

    const size_t pad_len = 500, N = data.size();
    kfr::univector<double> padded(N + 2 * pad_len), onda_lenta(N);
    std::copy(data.begin(), data.end(), padded.begin() + pad_len);
    const double left_edge_x2 = 2.0 * data[0];
    const double right_edge_x2 = 2.0 * data[N - 1];
    for (size_t i = 0; i < pad_len; i++)
    {
        padded[pad_len - 1 - i] = left_edge_x2 - data[i + 1];
        padded[N + pad_len + i] = right_edge_x2 - data[N - 2 - i];
    }
    kfr::filtfilt(padded, kfr::to_sos<double>(kfr::iir_lowpass(kfr::butterworth(4), fc, fs)));
    std::copy(padded.begin() + pad_len, padded.end() - pad_len, onda_lenta.begin());
    kfr::univector<double> onda_rapida = data - onda_lenta;

    std::ofstream out("results.csv");
    out << "t,original,onda_lenta,onda_rapida\n";
    for (size_t i = 0; i < data.size(); ++i)
    {
        out << (i / fs) << "," << data[i] << "," << onda_lenta[i] << "," << onda_rapida[i] << "\n";
    }

    return 0;
}