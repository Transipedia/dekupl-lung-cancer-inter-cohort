#include <iostream>
#include <string>
#include <unordered_map>
#include <set>
#include <fstream>
#include <cassert>
#include <algorithm>
#include <ctime>
#include <getopt.h>

#define MIN_ASSEMBLY_K 15

void complement(std::string &seq)
{
    auto lambda = [](const char c) {
        switch (c)
        {
        case 'A':
            return 'T';
        case 'G':
            return 'C';
        case 'C':
            return 'G';
        case 'T':
            return 'A';
        default:
            throw std::domain_error("Invalid nucleotide.");
        }
    };

    std::transform(seq.cbegin(), seq.cend(), seq.begin(), lambda);
}

void Contig_Anchor(std::unordered_map<std::string, std::set<std::string>> &kmer_index,
                   std::string &&contig_list_path,
                   const size_t k_length,
                   const bool stranded)
{
    std::ifstream cl_in(contig_list_path);
    std::string contig_name, line;
    size_t n_line = 0;

    if (!cl_in.is_open())
    {
        std::cerr << "ERROR: Contig list file " << contig_list_path << " does not exist, please use -h parameter for help..." << std::endl;
        exit(EXIT_FAILURE);
    }

    std::cerr << "Mapping contig list " << contig_list_path << "..." << std::endl;

    for (std::getline(cl_in, line); !cl_in.eof(); std::getline(cl_in, line))
    {
        if (line[0] == '>')
        {
            contig_name = line;
            contig_name.erase(0, 1);
            continue;
        }

        if (n_line > 0 && n_line % 5000 == 0)
        {
            std::cerr << "\t" << n_line << " contigs have been mapped." << std::endl;
        }

        if (line.size() < k_length)
        {
            continue;
        }

        for (size_t begin_pos = 0; begin_pos < line.size() - k_length + 1; begin_pos++)
        {
            auto kmer = line.substr(begin_pos, k_length);
            if (!stranded)
            {
                auto kmer_rc = kmer;
                reverse(kmer_rc.begin(), kmer_rc.end());
                complement(kmer_rc);
                if (kmer > kmer_rc)
                {
                    kmer = kmer_rc;
                }
            }

            auto ins_pair = kmer_index.insert({kmer, std::set<std::string>()});
            ins_pair.first->second.insert(contig_name);
        }
        n_line++;
    }
    cl_in.close();
}

void Print_KMer_Hash(const std::unordered_map<std::string, std::set<std::string>> &kmer_hash)
{
    for (const auto &elem : kmer_hash)
    {
        std::cerr << elem.first;
        for (const auto &elem2 : elem.second)
        {
            std::cerr << "\t" << elem2;
        }
        std::cerr << std::endl;
    }
}

void Make_Adjacency_List(std::unordered_map<std::string, int> &adjacency_list,
                         const std::unordered_map<std::string, std::set<std::string>> &kmer_hash)
{
    assert(adjacency_list.empty());

    for (const auto &elem : kmer_hash)
    {
        for (std::string contig_1 : elem.second)
        {
            for (std::string contig_2 : elem.second)
            {
                if (contig_1 >= contig_2)
                {
                    continue;
                }
                std::string contig_pair = contig_1 + "\t" + contig_2;
                auto iter_adjL = adjacency_list.find(contig_pair);
                if (iter_adjL == adjacency_list.cend())
                {
                    adjacency_list.insert({contig_pair, 1});
                }
                else
                {
                    iter_adjL->second++;
                }
            }
        }
    }
}

void Make_Adjacency_List2(std::unordered_map<std::string, int> &adjacency_list,
                          const std::unordered_map<std::string, std::set<std::string>> &kmer_hash_A,
                          const std::unordered_map<std::string, std::set<std::string>> &kmer_hash_B)
{
    assert(adjacency_list.empty());

    for (const auto &elem : kmer_hash_A)
    {
        std::string k_mer = elem.first;
        auto iter_B = kmer_hash_B.find(k_mer);
        if (iter_B == kmer_hash_B.cend())
        {
            continue;
        }
        for (std::string contig_A : elem.second)
        {
            for (std::string contig_B : iter_B->second)
            {
                std::string contig_pair = contig_A + "\t" + contig_B;
                auto iter_adjL = adjacency_list.find(contig_pair);
                if (iter_adjL == adjacency_list.cend())
                {
                    adjacency_list.insert({contig_pair, 1});
                }
                else
                {
                    iter_adjL->second++;
                }
            }
        }
    }
}

void Print_Adjacency_List(const std::unordered_map<std::string, int> &adjacency_list)
{
    std::cout << "contig_in_A\tcontig_in_B\tcount" << std::endl;

    for (const auto &elem : adjacency_list)
    {
        std::cout << elem.first << "\t" << elem.second << std::endl;
    }
}

int main(int argc, char *argv[])
{
    time_t t_begin = clock();
    std::string cohort_A_path, cohort_B_path;
    size_t k_length = 31;
    bool stranded = true;
    std::unordered_map<std::string, std::set<std::string>> kmer_hash;
    std::unordered_map<std::string, int> adjacency_list; // contig_name_1&congit_name_2, num //

    int c;
    while ((c = getopt(argc, argv, "nk:A:B:h")) >= 0)
    {
        switch (c)
        {
        case 'n':
            stranded = 0;
            break;
        case 'k':
            k_length = atoi(optarg);
            break;
        case 'A':
            cohort_A_path = optarg;
            break;
        case 'B':
            cohort_B_path = optarg;
            break;
        case 'h':
            std::cerr << std::endl;
            std::cerr << "Usage:\tconSeqGraph [-n] [-k <k_length>] -A cohort_list_A -B cohort_list_B" << std::endl;
            std::cerr << "Options:\t-k INT\t\tlength of k-mers (max_value: 32) [" << k_length << "]" << std::endl;
            std::cerr << "\t\t-n\t\tunstranded merging procedure" << std::endl;
            std::cerr << "\t\t-A String\tcohort A contig list path" << std::endl;
            std::cerr << "\t\t-B String\tcohort B contig list path" << std::endl;
            std::cerr << std::endl;
            exit(EXIT_FAILURE);
        }
    }

    std::cerr << "Stranded mode: " << (stranded ? "true" : "false") << std::endl;
    std::cerr << "K length: " << k_length << std::endl;
    std::cerr << "Cohort A path: " << cohort_A_path << std::endl;
    std::cerr << "Cohort B path: " << cohort_B_path << std::endl;

    Contig_Anchor(kmer_hash, std::move(cohort_A_path), k_length, stranded);
    Contig_Anchor(kmer_hash, std::move(cohort_B_path), k_length, stranded);
    Make_Adjacency_List(adjacency_list, kmer_hash);
    Print_Adjacency_List(adjacency_list);

    time_t t_end = clock();
    std::cerr << "==> Running time: " << static_cast<double>(t_end - t_begin) / CLOCKS_PER_SEC << " seconds." << std::endl;

    return EXIT_SUCCESS;
}
