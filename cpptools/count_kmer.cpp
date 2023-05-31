#include <iostream>
#include <string>
#include <fstream>
#include <map>
#include <unordered_map>
#include <cmath>
#include "cmdline.h"
#include "gzstream.h"
#include "ThreadPool/ThreadPool.h"
#include "util.h"

u_int64_t revComp(u_int64_t& x, int& sizeKmer)
{
    u_int64_t res = x;
    res = ((res >> 2 & 0x3333333333333333) | (res & 0x3333333333333333) << 2);
    res = ((res >> 4 & 0x0F0F0F0F0F0F0F0F) | (res & 0x0F0F0F0F0F0F0F0F) << 4);
    res = ((res >> 8 & 0x00FF00FF00FF00FF) | (res & 0x00FF00FF00FF00FF) << 8);
    res = ((res >> 16 & 0x0000FFFF0000FFFF) | (res & 0x0000FFFF0000FFFF) << 16);
    res = ((res >> 32 & 0x00000000FFFFFFFF) | (res & 0x00000000FFFFFFFF) << 32);
    res = res ^ 0xAAAAAAAAAAAAAAAA;
    return (res >> (2 * (32 - sizeKmer)));
}


std::pair<std::string, std::vector<double>> countKmer(
    int mlen, std::string reads_seq, std::string barcode, int k, int vector_size, int bin_size,
    std::unordered_map<unsigned long long, unsigned long>* kmer2frequency
)
{
    std::vector<double> hist(vector_size);

    if (barcode.empty() || reads_seq.size() <= mlen)
        return std::pair<std::string, std::vector<double>>("", hist);

    unsigned long total = 0;
    std::size_t len = 0;
    u_int64_t val = 0;
    unsigned long bases = (unsigned long)std::pow(2, k * 2) - 1;

    for (std::size_t i = 0; i < reads_seq.length(); ++i)
    {
        if (!(reads_seq[i] == 'A' || reads_seq[i] == 'C' || reads_seq[i] == 'G' || reads_seq[i] == 'T'))
        {
            val = 0;
            len = 0;
            continue;
        }
        val = (val << 2);
        val = val & bases;
        val += (reads_seq[i] >> 1 & 3);
        ++len;
        if (len == k)
        {
            --len;
            unsigned long long key = std::min(val, revComp(val, k));
            if (kmer2frequency->find(key) != kmer2frequency->end())
            {
                unsigned long count = kmer2frequency->at(key);
                int pos = count / bin_size;
                if (pos < vector_size)
                {
                    ++hist[pos];
                    ++total;
                }
            }
        }
    }

    // if (total)
    //     for (auto&& it : hist)
    //         it = it / total;
    // else
    //     for (auto&& it : hist)
    //         it = 1.0 / vector_size;

    return std::pair<std::string, std::vector<double>>(barcode, hist);
}

// count tnf
std::pair<std::string, std::vector<double>> countTnf(int mlen, std::string reads_seq, std::string barcode, int k, std::unordered_map<unsigned long long, unsigned long> kmer2frequency)
{
    std::vector<double> ret_fq;
    if (barcode.empty() || reads_seq.size() <= mlen)
        return std::pair<std::string, std::vector<double>>("", ret_fq);

    double total = 0;
    std::size_t len = 0;
    u_int64_t val = 0;
    unsigned long bases = (unsigned long)std::pow(2, k * 2) - 1;

    for (std::size_t i = 0; i < reads_seq.length(); ++i)
    {
        if (!(reads_seq[i] == 'A' || reads_seq[i] == 'C' || reads_seq[i] == 'G' || reads_seq[i] == 'T'))
        {
            val = 0;
            len = 0;
            continue;
        }
        val = (val << 2);
        val = val & bases;
        val += (reads_seq[i] >> 1 & 3);
        ++len;
        if (len == k)
        {
            --len;
            kmer2frequency[std::min(val, revComp(val, k))]++;
            total++;
        }
    }
    for (auto&& it : kmer2frequency)
        ret_fq.emplace_back(it.second);// / std::max(1.0, total));
    //ret_fq.emplace_back(it.second);

    return std::pair<std::string, std::vector<double>>(barcode, ret_fq);
}

std::pair<std::string, std::vector<double>> countTnfLongReads(int mlen, std::string reads_seq, std::string barcode, int k, std::unordered_map<unsigned long long, unsigned long> kmer2frequency)
{
    std::vector<double> ret_fq;
    if (barcode.empty() || reads_seq.size() <= mlen)
        return std::pair<std::string, std::vector<double>>("", ret_fq);

    double total = 0;
    std::size_t len = 0;
    u_int64_t val = 0;
    unsigned long bases = (unsigned long)std::pow(2, k * 2) - 1;

    for (std::size_t i = 0; i < reads_seq.length(); ++i)
    {
        if (!(reads_seq[i] == 'A' || reads_seq[i] == 'C' || reads_seq[i] == 'G' || reads_seq[i] == 'T'))
        {
            val = 0;
            len = 0;
            continue;
        }
        val = (val << 2);
        val = val & bases;
        val += (reads_seq[i] >> 1 & 3);
        ++len;
        if (len == k)
        {
            --len;
            kmer2frequency[std::min(val, revComp(val, k))]++;
            total++;
        }
    }
    for (auto&& it : kmer2frequency)
        ret_fq.emplace_back(it.second);// / std::max(1.0, total));
    //ret_fq.emplace_back(it.second);

    return std::pair<std::string, std::vector<double>>(barcode, ret_fq);
}


std::pair<std::string, std::vector<double>> countKmerLongReads(
    int mlen, std::string reads_seq, std::string read_name, int k, int vector_size, int bin_size,
    std::unordered_map<unsigned long long, unsigned long>* kmer2frequency
)
{
    std::vector<double> hist(vector_size);

    if (read_name.empty() || reads_seq.size() <= mlen)
        return std::pair<std::string, std::vector<double>>("", hist);
    unsigned long total = 0;
    std::size_t len = 0;
    u_int64_t val = 0;
    unsigned long bases = (unsigned long)std::pow(2, k * 2) - 1;

    for (std::size_t i = 0; i < reads_seq.length(); ++i)
    {
        if (!(reads_seq[i] == 'A' || reads_seq[i] == 'C' || reads_seq[i] == 'G' || reads_seq[i] == 'T'))
        {
            val = 0;
            len = 0;
            continue;
        }
        val = (val << 2);
        val = val & bases;
        val += (reads_seq[i] >> 1 & 3);
        ++len;
        if (len == k)
        {
            --len;
            unsigned long long key = std::min(val, revComp(val, k));
            if (kmer2frequency->find(key) != kmer2frequency->end())
            {
                unsigned long count = kmer2frequency->at(key);
                int pos = count / bin_size;
                if (pos < vector_size)
                {
                    ++hist[pos];
                    ++total;
                }
            }
        }
    }

    // if (total)
    //     for (auto&& it : hist)
    //         it = it / total;
    // else
    //     for (auto&& it : hist)
    //         it = 1.0 / vector_size;

    return std::pair<std::string, std::vector<double>>(read_name, hist);
}

int getKmer2Freq_kmer(const std::string global, std::unordered_map<unsigned long long, unsigned long>& kmer2frequency, int kmer)
{
    std::string line;
    std::fstream globalf(global.c_str(), std::fstream::in);
    unsigned long bases = (unsigned long)std::pow(2, kmer * 2) - 1;

    while (getline(globalf, line))
    {
        std::size_t pos = line.find_first_of('\t');
        std::string it = line.substr(0, pos);
        unsigned long freq = std::stol(line.substr(pos + 1));
        std::size_t len = 0;
        u_int64_t val = 0;
        for (std::size_t i = 0; i < it.length(); ++i)
        {
            if (!(it[i] == 'A' || it[i] == 'C' || it[i] == 'G' || it[i] == 'T'))
            {
                val = 0;
                len = 0;
                continue;
            }
            val = (val << 2);
            val = val & bases;
            val += (it[i] >> 1 & 3);
            ++len;
            if (len == kmer)
            {
                --len;
                kmer2frequency[std::min(val, revComp(val, kmer))] = freq;
            }
        }
    }
    globalf.close();
    return 0;
}

std::vector<std::string> generateAllTNFs(int n)
{
    std::vector<std::string> ret;
    if (n == 1)
    {
        ret.emplace_back("A");
        ret.emplace_back("T");
        ret.emplace_back("C");
        ret.emplace_back("G");
    }
    else
    {
        std::vector<std::string> ret_last = std::move(generateAllTNFs(n - 1));
        for (auto&& kmer : ret_last)
        {
            ret.emplace_back("A" + kmer);
            ret.emplace_back("T" + kmer);
            ret.emplace_back("C" + kmer);
            ret.emplace_back("G" + kmer);
        }
    }
    return ret;
}

int getTnf2Frequency(std::unordered_map<unsigned long long, unsigned long>& kmer2frequency, int kmer)
{
    std::vector<std::string> all_kmers = std::move(generateAllTNFs(kmer));
    unsigned long bases = (unsigned long)std::pow(2, kmer * 2) - 1;

    for (auto&& it : all_kmers)
    {
        std::size_t len = 0;
        u_int64_t val = 0;
        for (std::size_t i = 0; i < it.length(); ++i)
        {
            if (!(it[i] == 'A' || it[i] == 'C' || it[i] == 'G' || it[i] == 'T'))
            {
                val = 0;
                len = 0;
                continue;
            }
            val = (val << 2);
            val = val & bases;
            val += (it[i] >> 1 & 3);
            len++;
            if (len == kmer)
            {
                len--;
                kmer2frequency[std::min(val, revComp(val, kmer))] = 0;
            }
        }
    }
    return 0;
}



int main(int argc, char* argv[])
{
    cmdline::parser argParser;
    argParser.add<std::string>("reads1", '1', "Path to reads1 (gzipped and barcode sorted).", false);
    argParser.add<std::string>("reads2", '2', "Path to reads2 (gzipped and barcode sorted).", false);
    argParser.add<std::string>("global", 'g', "Path to global frequency by jellyfish.", false);
    argParser.add<std::string>("output", 'o', "Output (gzipped).", true);
    argParser.add<int>("kmer", 'k', "Length of k.", false, 15);
    argParser.add<int>("len", 'l', "Min length of barcode.", false, 1000);
    argParser.add<int>("thread", 't', "Number of thread to use.", false, 16);
    argParser.add<int>("vector", 'v', "Vector size of frequency.", false, 400);
    argParser.add<int>("window", 'w', "Window size of frequency.", false, 10);
    argParser.add<std::string>("interleaved", 'i', "Path to interleaved reads. If specified, --reads1 and --reads2 are not needed.", false);
    argParser.add<std::string>("long_reads", 'r', "Path to long reads. If specified, --reads1 and --reads2 are not needed.", false);
    
    argParser.parse_check(argc, argv);

    std::string reads1 = argParser.get<std::string>("reads1");
    std::string reads2 = argParser.get<std::string>("reads2");
    std::string interleaved = argParser.get<std::string>("interleaved");
    std::string long_reads = argParser.get<std::string>("long_reads");

    std::string global = argParser.get<std::string>("global");
    std::string output = argParser.get<std::string>("output");
    int kmer = argParser.get<int>("kmer");
    int mlen = argParser.get<int>("len");
    int thread = argParser.get<int>("thread");
    int vector_size = argParser.get<int>("vector");
    int bin_size = argParser.get<int>("window");

    ThreadPool thread_pool(thread);
    std::vector<std::future<std::pair<std::string, std::vector<double>>>> results;

    std::unordered_map<unsigned long long, unsigned long> kmer2frequency;

    if (kmer == 15) {
        getKmer2Freq_kmer(global, kmer2frequency, kmer);
    }
    else if (kmer == 4 || kmer == 3) {
        getTnf2Frequency(kmer2frequency, kmer);
    }
    
    ogzstream outputf(output.c_str());
    if (!reads1.empty() && !reads2.empty())
    {
        igzstream reads1f(reads1.c_str()), reads2f(reads2.c_str());
        unsigned long long cnt_line = 0, cnt_unpaired = 0;
        std::string line1, line2, last_barcode, reads_seq;
        std::pair<std::string, std::string> p1, p2;

        while (getline(reads1f, line1))
        {
            getline(reads2f, line2);
            switch (++cnt_line % 4)
            {
            case 1:
                p1 = std::move(getBarcode(line1)), p2 = std::move(getBarcode(line2));
                break;
            case 2:
                if (p1.first.compare(p2.first) || p1.second.compare(p2.second))
                    ++cnt_unpaired;
                else
                {
                    reads_seq += (line1 + "N" + line2 + "N");
                    if (p1.second.compare(last_barcode))
                    {
                        if (results.size() >= 10000)
                        {
                            for (auto&& result : results)
                            {
                                std::pair<std::string, std::vector<double>> ret = std::move(result.get());
                                if (ret.first.empty())
                                    continue;
                                outputf << ret.first;
                                for (auto&& fq : ret.second)
                                    outputf << ',' << fq;
                                outputf << std::endl;
                            }
                            results.clear();
                        }
                        results.emplace_back(thread_pool.enqueue(countKmer, mlen, reads_seq, last_barcode, kmer, vector_size, bin_size, &kmer2frequency));
                        last_barcode = std::move(p1.second);
                        reads_seq.clear();
                    }
                }
                break;
            case 0:
                if (cnt_line % 4000000 == 0)
                    std::cout << "kmer Processed " << cnt_line / 4 << " read pairs. " << cnt_unpaired << " unpaired." << std::endl;
                break;
            default:
                break;
            }
        }
        results.emplace_back(thread_pool.enqueue(countKmer, mlen, reads_seq, last_barcode, kmer, vector_size, bin_size, &kmer2frequency));
        reads1f.close();
        reads2f.close();
    }
    else if (!interleaved.empty())
    {
        std::cout << "interleaved file: "<<interleaved.c_str()<<std::endl;
        igzstream readsf(interleaved.c_str());
        unsigned long long cnt_line = 0, cnt_unpaired = 0;
        std::string line, last_barcode, barcode, reads_seq;
        std::size_t pos1, pos2;

        while (getline(readsf, line))
        {
            switch (++cnt_line % 8)
            {
            case 1:
                // pos1 = line.find_first_of('\t');
                // barcode.clear();
                // if (pos1 != std::string::npos)
                // {
                //     pos2 = line.find_first_of('-', pos1 + 6);
                //     barcode = line.substr(pos1 + 6, pos2 - pos1 - 6);
                // }
                barcode = getBarcode(line).second;
                break;
            case 2:
                reads_seq += (line + "N\n");
                break;
            case 6:
                reads_seq += (line + "N\n");
                if (barcode.compare(last_barcode))
                {
                    if (results.size() >= 10000)
                    {
                        for (auto&& result : results)
                        {
                            std::pair<std::string, std::vector<double>> ret = std::move(result.get());
                            if (ret.first.empty())
                                continue;
                            outputf << ret.first;
                            for (auto&& fq : ret.second)
                                outputf << ',' << fq;
                            outputf << std::endl;
                        }
                        results.clear();
                    }
                    if (kmer == 15) {
                        results.emplace_back(thread_pool.enqueue(countKmer, mlen, reads_seq, last_barcode, kmer, vector_size, bin_size, &kmer2frequency));
                    } else if (kmer == 4) {
                        results.emplace_back(thread_pool.enqueue(countTnf, mlen, reads_seq, last_barcode, kmer, kmer2frequency));
                    }
                    last_barcode = std::move(barcode);
                    reads_seq.clear();
                }
                break;
            case 0:
                if (cnt_line % 8000000 == 0)
                    std::cout << "kmer Processed " << cnt_line / 8 << " read pairs." << std::endl;
                break;
            default:
                break;
            }
        }
        if (kmer == 15) {
            results.emplace_back(thread_pool.enqueue(countKmer, mlen, reads_seq, last_barcode, kmer, vector_size, bin_size, &kmer2frequency));
        } else if (kmer == 4) {
            results.emplace_back(thread_pool.enqueue(countTnf, mlen, reads_seq, last_barcode, kmer, kmer2frequency));
        }
        readsf.close();
    }
    else if (!long_reads.empty())
    {
        igzstream longreadsf(long_reads.c_str());
        unsigned long long cnt_line = 0, cnt_unpaired = 0;
        std::string line, read_name, reads_seq;
        // std::pair<std::string, std::string> read_line;

        while (getline(longreadsf, line))
        {
            switch (++cnt_line % 4)
            {
            case 1:
                read_name = std::move(getReadname(line));
                break;
            case 2:
                reads_seq = line;
                if (results.size() >= 10000)
                {
                    for (auto&& result : results)
                    {
                        std::pair<std::string, std::vector<double>> ret = std::move(result.get());
                        if (ret.first.empty())
                            continue;
                        outputf << ret.first;
                        for (auto&& fq : ret.second)
                            outputf << ',' << fq;
                        outputf << std::endl;
                    }
                    results.clear();
                }
                if (kmer == 15) {
                    results.emplace_back(thread_pool.enqueue(countKmerLongReads, mlen, reads_seq, read_name, kmer, vector_size, bin_size, &kmer2frequency));
                } else if (kmer == 4 || kmer == 3) {
                    results.emplace_back(thread_pool.enqueue(countTnfLongReads, mlen, reads_seq, read_name, kmer, kmer2frequency));
                }
                reads_seq.clear();
                
                break;
            case 0:
                if (cnt_line % 8000000 == 0)
                    std::cout << "kmer Processed " << cnt_line / 4 << " read pairs. " << cnt_unpaired << " unpaired." << std::endl;
                break;
            default:
                break;
            }
        }
    }
    for (auto&& result : results)
    {
        std::pair<std::string, std::vector<double>> ret = std::move(result.get());
        if (ret.first.empty())
            continue;
        outputf << ret.first;
        for (auto&& fq : ret.second)
            outputf << ',' << fq;
        outputf << std::endl;
    }
    outputf.close();
    return 0;
}