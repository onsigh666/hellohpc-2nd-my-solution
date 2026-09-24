#include "main.hpp"

#include <time.h>

#include <fstream>
#include <iostream>

#include "timer.hpp"

using namespace std;

const uint32 MD5IV[] = {0x67452301, 0xefcdab89, 0x98badcfe, 0x10325476};

unsigned load_block(istream& i, uint32 block[]);
void save_block(ostream& o, const uint32 block[]);
// 返回 false 表示第一块搜索被 g_stop 叫停，输出参数全部作废。
bool find_collision(const uint32 IV[], uint32 msg1block0[], uint32 msg1block1[],
                    uint32 msg2block0[], uint32 msg2block1[],
                    bool verbose = false);

#if 0

// example trivial version with md5 initial value
int main()
{
	seed32_1 = uint32(time(NULL));
	seed32_2 = 0x12345678;

	uint32 IV[4] = { MD5IV[0], MD5IV[1], MD5IV[2], MD5IV[3] };
	uint32 msg1block0[16];
	uint32 msg1block1[16];
	uint32 msg2block0[16];
	uint32 msg2block1[16];
	find_collision(IV, msg1block0, msg1block1, msg2block0, msg2block1, true);

	ofstream ofs1("msg1", ios::binary);
	save_block(ofs1, msg1block0);
	save_block(ofs1, msg1block1);
	ofstream ofs2("msg2", ios::binary);
	save_block(ofs2, msg2block0);
	save_block(ofs2, msg2block1);
}

#else

#include <ctype.h>
#include <errno.h>
#include <sched.h>
#include <stdint.h>
#include <stdlib.h>
#include <sys/stat.h>

#include <array>
#include <atomic>
#include <memory>
#include <mutex>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <utility>
#include <vector>

typedef ::uint64_t uint64;

namespace {

const char* const options_help =
    "Allowed options:\n"
    "  -h [ --help ]           Show options.\n"
    "  -q [ --quiet ]          Be less verbose.\n"
    "  -i [ --ihv ] arg        Use specified initial value. Default is MD5 "
    "initial \n"
    "                          value.\n"
    "  --seed1 arg             Specify a SEED1 value.  Default is based on "
    "time.\n"
    "  --seed2 arg             Specify a SEED2 value.  Default is a constant "
    "value.\n"
    "  -p [ --prefixfile ] arg Calculate initial value using given prefixfile. "
    "Also \n"
    "                          copies data to output files.\n"
    "  -o [ --out ] arg        Set output filenames. This must be the last "
    "option \n"
    "                          and exactly 2 filenames must be specified. \n"
    "                          Default: -o msg1.bin msg2.bin\n";

bool file_exists(const std::string& path) {
  struct stat st;
  return ::stat(path.c_str(), &st) == 0;
}

std::string uint_to_string(unsigned value) {
  std::ostringstream oss;
  oss << value;
  return oss.str();
}

uint32 parse_uint32(const std::string& name, const std::string& value) {
  // An optional sign followed by digits only; no surrounding whitespace.
  // Out-of-range values are rejected, while an in-range negative value wraps.
  std::string::size_type pos = 0;
  bool negative = false;
  if (pos < value.size() && (value[pos] == '+' || value[pos] == '-'))
    negative = (value[pos++] == '-');

  bool ok = pos < value.size();
  for (std::string::size_type i = pos; ok && i < value.size(); ++i)
    if (!isdigit((unsigned char)value[i])) ok = false;

  uint64 magnitude = 0;
  if (ok) {
    errno = 0;
    magnitude = strtoull(value.c_str() + pos, 0, 10);
    if (errno == ERANGE || magnitude > 0xFFFFFFFFull) ok = false;
  }

  if (!ok)
    throw std::runtime_error("the argument ('" + value + "') for option '--" +
                             name + "' is invalid");

  uint32 result = uint32(magnitude);
  return negative ? uint32(0) - result : result;
}

struct options {
  bool help;
  bool quiet;
  bool has_ihv;
  bool has_prefixfile;
  bool has_out;
  bool testmd5iv;
  bool testrndiv;
  bool testreciv;
  bool testall;
  bool has_tasks;
  std::string tasks;
  std::string ihv;
  std::string prefixfile;
  std::vector<std::string> out;

  options()
      : help(false),
        quiet(false),
        has_ihv(false),
        has_prefixfile(false),
        has_out(false),
        testmd5iv(false),
        testrndiv(false),
        testreciv(false),
        has_tasks(false),
        testall(false) {}
};

// Returns the argument for an option that requires one, either the remainder of
// the current token or the next argument.
std::string take_value(const std::string& name, const std::string& rest,
                       int argc, char** argv, int& i) {
  if (!rest.empty()) return rest;
  if (i + 1 >= argc)
    throw std::runtime_error("the required argument for option '--" + name +
                             "' is missing");
  return argv[++i];
}

void parse_command_line(int argc, char** argv, options& opt, uint32& seed32_1,
                        uint32& seed32_2) {
  bool have_positional = false;

  for (int i = 1; i < argc; ++i) {
    std::string arg = argv[i];

    if (arg.size() > 2 && arg.compare(0, 2, "--") == 0) {
      std::string name = arg.substr(2);
      std::string value;
      bool has_value = false;
      std::string::size_type eq = name.find('=');
      if (eq != std::string::npos) {
        value = name.substr(eq + 1);
        name = name.substr(0, eq);
        has_value = true;
      }

      if (name == "help")
        opt.help = true;
      else if (name == "quiet")
        opt.quiet = true;
      else if (name == "testmd5iv")
        opt.testmd5iv = true;
      else if (name == "testrndiv")
        opt.testrndiv = true;
      else if (name == "testreciv")
        opt.testreciv = true;
      else if (name == "testall")
        opt.testall = true;
      else if (name == "ihv") {
        opt.ihv = has_value ? value : take_value(name, "", argc, argv, i);
        opt.has_ihv = true;
      } else if (name == "prefixfile") {
        opt.prefixfile =
            has_value ? value : take_value(name, "", argc, argv, i);
        opt.has_prefixfile = true;
      } else if (name == "seed1")
        seed32_1 = parse_uint32(
            name, has_value ? value : take_value(name, "", argc, argv, i));
      else if (name == "seed2")
        seed32_2 = parse_uint32(
            name, has_value ? value : take_value(name, "", argc, argv, i));
      else if (name == "out") {
        // Multitoken: consume every following argument that is not an option.
        if (has_value)
          opt.out.push_back(value);
        else
          opt.out.push_back(take_value(name, "", argc, argv, i));
        while (i + 1 < argc && argv[i + 1][0] != '-')
          opt.out.push_back(argv[++i]);
        opt.has_out = true;
      } else if (name == "tasks") {
        opt.tasks = has_value ? value : take_value(name, "", argc, argv, i);
        opt.has_tasks = true;
      } else
        throw std::runtime_error("unrecognised option '" + arg + "'");
    } else if (arg.size() > 1 && arg[0] == '-') {
      // Short options may be grouped; those taking an argument consume the
      // rest of the token, or the next argument if the token ends there.
      for (std::string::size_type c = 1; c < arg.size(); ++c) {
        char o = arg[c];
        std::string rest = arg.substr(c + 1);
        if (o == 'h')
          opt.help = true;
        else if (o == 'q')
          opt.quiet = true;
        else if (o == 'i') {
          opt.ihv = take_value("ihv", rest, argc, argv, i);
          opt.has_ihv = true;
          break;
        } else if (o == 'p') {
          opt.prefixfile = take_value("prefixfile", rest, argc, argv, i);
          opt.has_prefixfile = true;
          break;
        } else if (o == 'o') {
          opt.out.push_back(take_value("out", rest, argc, argv, i));
          while (i + 1 < argc && argv[i + 1][0] != '-')
            opt.out.push_back(argv[++i]);
          opt.has_out = true;
          break;
        } else
          throw std::runtime_error(std::string("unrecognised option '-") + o +
                                   "'");
      }
    } else {
      // Sole positional argument is the prefixfile.
      if (have_positional || opt.has_prefixfile)
        throw std::runtime_error(
            "too many positional options have been "
            "specified on the command line");
      opt.prefixfile = arg;
      opt.has_prefixfile = true;
      have_positional = true;
    }
  }
}

}  // namespace

/* ---------- 批量模式：一个进程内跑完 tasks.txt ---------- */

// 把 64 位标签打散成 32 位种子，并强制最低位为 1。
// 为什么必须非 0：xorshift 的 0 是"吸收态"——种子为 0 时永远只产生 0，
// 搜索会退化成死循环。mix 用的是 splitmix64 的尾部，雪崩性够好。
static uint32 mix32(uint64 x) {
  x += 0x9E3779B97F4A7C15ull;
  x = (x ^ (x >> 30)) * 0xBF58476D1CE4E5B9ull;
  x = (x ^ (x >> 27)) * 0x94D049BB133111EBull;
  return uint32(x ^ (x >> 31)) | 1u;
}

// 读前缀、链式压缩求出 IV。不写任何文件 —— 现在有 P 个搜索在抢同一个
// 任务，只有赢家才写输出，所以算 IV 和写输出必须分开。
static bool prefix_iv(const std::string& infile, uint32 IV[4]) {
  IV[0] = MD5IV[0];
  IV[1] = MD5IV[1];
  IV[2] = MD5IV[2];
  IV[3] = MD5IV[3];

  std::ifstream ifs(infile.c_str(), std::ios::binary);
  if (!ifs) return false;

  uint32 block[16];
  while (true) {
    unsigned len = load_block(ifs, block);
    if (!len) break;
    md5_compress(IV, block);
  }
  return true;
}

// 只有赢家调用：把前缀（补零到 64 字节整块）和 4 个碰撞块写进两个文件。
// load_block 会把不足 64 字节的最后一块补零，save_block 又把补的零一起
// 写出去 —— 所以输出文件里前缀被补到整块，碰撞块才对齐得上。
static bool write_outputs(const std::string& infile, const std::string& out1,
                          const std::string& out2, const uint32* m1b0,
                          const uint32* m1b1, const uint32* m2b0,
                          const uint32* m2b1) {
  std::ifstream ifs(infile.c_str(), std::ios::binary);
  if (!ifs) return false;
  std::ofstream ofs1(out1.c_str(), std::ios::binary);
  if (!ofs1) return false;
  std::ofstream ofs2(out2.c_str(), std::ios::binary);
  if (!ofs2) return false;

  uint32 block[16];
  while (true) {
    unsigned len = load_block(ifs, block);
    if (!len) break;
    save_block(ofs1, block);
    save_block(ofs2, block);
  }

  save_block(ofs1, m1b0);
  save_block(ofs1, m1b1);
  save_block(ofs2, m2b0);
  save_block(ofs2, m2b1);

  ofs1.flush();
  ofs2.flush();
  return ofs1.good() && ofs2.good();
}

// 本进程实际能用几个核。顺序很重要：
//   1. OMP_NUM_THREADS —— 评测脚本/集群 sbatch 会设，优先级最高
//   2. CPU 亲和掩码   —— sched_getaffinity 会尊重 cgroup/cpuset 限制
//   3. 兜底            —— hardware_concurrency 只看整机，作业被限核时会超订
static unsigned usable_cores() {
  const char* e = getenv("OMP_NUM_THREADS");
  if (e && *e) {
    int n = atoi(e);
    if (n > 0) return unsigned(n);
  }
  cpu_set_t set;
  if (sched_getaffinity(0, sizeof(set), &set) == 0) {
    int n = CPU_COUNT(&set);
    if (n > 0) return unsigned(n);
  }
  unsigned n = std::thread::hardware_concurrency();
  return n ? n : 1u;
}

int run_tasks(const char* taskfn) {
  std::ifstream tf(taskfn);
  if (!tf) {
    std::cerr << "cannot read task file: " << taskfn << "\n";
    return 1;
  }

  std::vector<std::array<std::string, 3> > tasks;
  std::string a, b, c;
  while (tf >> a >> b >> c) {
    std::array<std::string, 3> t = {a, b, c};
    tasks.push_back(t);
  }
  if (tasks.empty()) {
    std::cerr << "no tasks\n";
    return 1;
  }

  const size_t N = tasks.size();
  const unsigned P = usable_cores();
  const size_t K = P;       // 每个任务同时开几个搜索 —— 就按核数来
  const size_t NM = N * K;  // "搜索"总数：任务数 × K

  // std::atomic 不可拷贝，所以用 unique_ptr 数组
  std::unique_ptr<std::atomic<bool>[]> solved(new std::atomic<bool>[N]);
  for (size_t i = 0; i < N; ++i) solved[i].store(false);

  std::cerr << "tasks=" << N << " threads=" << P << " missions=" << NM << "\n";

  std::atomic<size_t> next(0);  // 下一个该领的搜索编号
  std::atomic<int> failures(0);
  std::mutex cout_mtx;

  std::vector<std::thread> pool;
  for (unsigned t = 0; t < P; ++t)
    pool.push_back(std::thread([&, t]() {
      size_t m;
      while ((m = next.fetch_add(1)) < NM)  // 动态领搜索
      {
        const size_t i = m / K;  // 第 m 个搜索属于第几个任务
        if (solved[i].load(std::memory_order_relaxed))
          continue;  // 这个任务已经被人搞定了

        uint64 tag = uint64(m) * 0x9E3779B97F4A7C15ull + t;
        uint32 s1 = mix32(tag), s2 = mix32(tag + 0x1234567ull);

        uint32 IV[4];
        if (!prefix_iv(tasks[i][0], IV)) {
          ++failures;
          continue;
        }

        uint32 m1b0[16], m1b1[16], m2b0[16], m2b1[16];
        g_stop = &solved[i];  // 装好叫停开关
        seed32_1 = s1;
        seed32_2 = s2;
        bool found = find_collision(IV, m1b0, m1b1, m2b0, m2b1);
        g_stop = nullptr;

        if (!found) continue;  // 被叫停，我的结果作废

        // exchange 返回 true 表示之前已经是 true —— 别人先赢了
        if (solved[i].exchange(true)) continue;

        bool ok = write_outputs(tasks[i][0], tasks[i][1], tasks[i][2], m1b0,
                                m1b1, m2b0, m2b1);
        if (!ok) ++failures;
        std::lock_guard<std::mutex> lk(cout_mtx);
        std::cout << (ok ? "ok " : "FAIL ") << i << " " << tasks[i][0] << "\n"
                  << std::flush;
      }
    }));
  for (unsigned t = 0; t < P; ++t) pool[t].join();

  std::cerr << "done, failures=" << failures.load() << "\n";
  return failures.load() ? 1 : 0;
}

void test_md5iv(bool single = false);
void test_rndiv(bool single = false);
void test_reciv(bool single = false);
void test_all();
int run_tasks(const char* taskfn);

int main(int argc, char** argv) {
  seed32_1 = uint32(time(NULL));
  seed32_2 = 0x12345678;

  uint32 IV[4] = {MD5IV[0], MD5IV[1], MD5IV[2], MD5IV[3]};

  string outfn1 = "msg1.bin";
  string outfn2 = "msg2.bin";
  string ihv;
  string prefixfn;
  bool verbose = true;

  cout << "MD5 collision generator v1.5\n"
          "by Marc Stevens (http://www.win.tue.nl/hashclash/)\n"
       << endl;

  try {
    hashclash::timer runtime(true);

    options opt;
    parse_command_line(argc, argv, opt, seed32_1, seed32_2);
    ihv = opt.ihv;
    prefixfn = opt.prefixfile;

    if (opt.quiet) verbose = false;

    if (opt.has_tasks) return run_tasks(opt.tasks.c_str());

    if (opt.help || argc == 1) {
      cout << options_help << endl;
      return 1;
    }

    if (opt.testmd5iv) test_md5iv();
    if (opt.testrndiv) test_rndiv();
    if (opt.testreciv) test_reciv();
    if (opt.testall) test_all();

    if (opt.has_prefixfile) {
      unsigned l = prefixfn.size();
      if (l >= 4 && prefixfn[l - 4] == '.' && prefixfn[l - 3] != '.' &&
          prefixfn[l - 2] != '.' && prefixfn[l - 1] != '.') {
        outfn1 = prefixfn.substr(0, l - 4) + "_msg1" + prefixfn.substr(l - 4);
        outfn2 = prefixfn.substr(0, l - 4) + "_msg2" + prefixfn.substr(l - 4);
        unsigned i = 1;
        while (file_exists(outfn1) || file_exists(outfn2)) {
          outfn1 = prefixfn.substr(0, l - 4) + "_msg1_" + uint_to_string(i) +
                   prefixfn.substr(l - 4);
          outfn2 = prefixfn.substr(0, l - 4) + "_msg2_" + uint_to_string(i) +
                   prefixfn.substr(l - 4);
          ++i;
        }
      }
    }

    if (opt.has_out) {
      vector<string>& outfns = opt.out;
      if (outfns.size() != 2) {
        cerr << "Error: exactly two output filenames should be specified."
             << endl;
        return 1;
      }
      outfn1 = outfns[0];
      outfn2 = outfns[1];
    }

    if (verbose)
      cout << "Using output filenames: '" << outfn1 << "' and '" << outfn2
           << "'" << endl;
    ofstream ofs1(outfn1.c_str(), ios::binary);
    if (!ofs1) {
      cerr << "Error: cannot open outfile: '" << outfn1 << "'" << endl;
      return 1;
    }
    ofstream ofs2(outfn2.c_str(), ios::binary);
    if (!ofs2) {
      cerr << "Error: cannot open outfile: '" << outfn2 << "'" << endl;
      return 1;
    }

    if (opt.has_prefixfile) {
      if (verbose) cout << "Using prefixfile: '" << prefixfn << "'" << endl;
      ifstream ifs(prefixfn.c_str(), ios::binary);
      if (!ifs) {
        cerr << "Error: cannot open inputfile: '" << prefixfn << "'" << endl;
        return 1;
      }
      uint32 block[16];
      while (true) {
        unsigned len = load_block(ifs, block);
        if (len) {
          save_block(ofs1, block);
          save_block(ofs2, block);
          md5_compress(IV, block);
        } else
          break;
      }
    } else {
      if (!opt.has_ihv) ihv = "0123456789abcdeffedcba9876543210";
      if (ihv.size() != 32) {
        cerr << "Error: an initial value must be specified as a hash value of "
                "32 hexadecimal characters."
             << endl;
        return 1;
      } else {
        uint32 c;
        for (unsigned i = 0; i < 4; ++i) {
          IV[i] = 0;
          for (unsigned b = 0; b < 4; ++b) {
            stringstream ss;
            ss << ihv.substr(i * 8 + b * 2, 2);
            ss >> hex >> c;
            IV[i] += c << (b * 8);
          }
        }
      }
    }
    if (verbose) {
      cout << "Using initial value: " << hex;
      unsigned oldwidth = cout.width(2);
      char oldfill = cout.fill('0');

      for (unsigned i = 0; i < 4; ++i) {
        for (unsigned b = 0; b < 4; ++b) {
          cout.width(2);
          cout.fill('0');
          cout << ((IV[i] >> (b * 8)) & 0xFF);
        }
      }
      cout.width(oldwidth);
      cout.fill(oldfill);
      cout << dec << endl;
    }

    if (verbose) cout << endl;

    uint32 msg1block0[16];
    uint32 msg1block1[16];
    uint32 msg2block0[16];
    uint32 msg2block1[16];
    find_collision(IV, msg1block0, msg1block1, msg2block0, msg2block1, true);

    save_block(ofs1, msg1block0);
    save_block(ofs1, msg1block1);
    save_block(ofs2, msg2block0);
    save_block(ofs2, msg2block1);
    if (verbose) cout << "Running time: " << runtime.time() << " s" << endl;
    return 0;
  } catch (exception& e) {
    cerr << "\nException caught:\n" << e.what() << endl;
    return 1;
  } catch (...) {
    cerr << "\nUnknown exception caught!" << endl;
    return 1;
  }
}

void test_md5iv(bool single) {
  uint32 IV[4] = {MD5IV[0], MD5IV[1], MD5IV[2], MD5IV[3]};
  uint32 msg1block0[16];
  uint32 msg1block1[16];
  uint32 msg2block0[16];
  uint32 msg2block1[16];

  hashclash::timer runtime(true);
  while (true) {
    runtime.start();
    find_collision(IV, msg1block0, msg1block1, msg2block0, msg2block1);
    double time = runtime.time();
    cout << endl << "Running time: " << time << " s" << endl;
    ofstream of_timings("timings_md5iv.txt", ios::app);
    of_timings << time << endl;
    if (single) return;
  }
}

void test_rndiv(bool single) {
  uint32 IV[4];
  uint32 msg1block0[16];
  uint32 msg1block1[16];
  uint32 msg2block0[16];
  uint32 msg2block1[16];

  hashclash::timer runtime(true);
  while (true) {
    runtime.start();
    IV[0] = xrng64();
    IV[1] = xrng64();
    IV[2] = xrng64();
    IV[3] = xrng64();
    find_collision(IV, msg1block0, msg1block1, msg2block0, msg2block1);
    double time = runtime.time();
    cout << endl << "Running time: " << time << " s" << endl;
    ofstream of_timings("timings_rndiv.txt", ios::app);
    of_timings << time << endl;
    if (single) return;
  }
}

void test_reciv(bool single) {
  uint32 IV[4];
  uint32 msg1block0[16];
  uint32 msg1block1[16];
  uint32 msg2block0[16];
  uint32 msg2block1[16];

  hashclash::timer runtime(true);
  while (true) {
    runtime.start();
    IV[0] = xrng64();
    IV[1] = xrng64();
    IV[2] = xrng64();
    IV[3] = xrng64();
    IV[2] |= 1 << 25;
    IV[2] ^= ((IV[2] & (1 << 24)) << 1);
    IV[3] &= ~(1 << 25);
    IV[3] ^= ((IV[3] & (1 << 24)) << 1);

    find_collision(IV, msg1block0, msg1block1, msg2block0, msg2block1);
    double time = runtime.time();
    cout << endl << "Running time: " << time << " s" << endl;
    ofstream of_timings("timings_reciv.txt", ios::app);
    of_timings << time << endl;
    if (single) return;
  }
}

void test_all() {
  while (true) {
    test_md5iv(true);
    test_reciv(true);
    test_rndiv(true);
  }
}

#endif

unsigned load_block(istream& i, uint32 block[]) {
  unsigned len = 0;
  char uc;
  for (unsigned k = 0; k < 16; ++k) {
    block[k] = 0;
    for (unsigned c = 0; c < 4; ++c) {
      i.get(uc);
      if (i)
        ++len;
      else
        uc = 0;
      block[k] += uint32((unsigned char)(uc)) << (c * 8);
    }
  }
  return len;
}

void save_block(ostream& o, const uint32 block[]) {
  for (unsigned k = 0; k < 16; ++k)
    for (unsigned c = 0; c < 4; ++c)
      o << (unsigned char)((block[k] >> (c * 8)) & 0xFF);
}

bool find_collision(const uint32 IV[], uint32 msg1block0[], uint32 msg1block1[],
                    uint32 msg2block0[], uint32 msg2block1[], bool verbose) {
  if (verbose) cout << "Generating first block: " << flush;
  if (!find_block0(msg1block0, IV)) return false;  // 被叫停，后面全部作废

  uint32 IHV[4] = {IV[0], IV[1], IV[2], IV[3]};
  md5_compress(IHV, msg1block0);

  if (verbose) cout << endl << "Generating second block: " << flush;
  find_block1(msg1block1, IHV);

  for (int t = 0; t < 16; ++t) {
    msg2block0[t] = msg1block0[t];
    msg2block1[t] = msg1block1[t];
  }
  msg2block0[4] += 1 << 31;
  msg2block0[11] += 1 << 15;
  msg2block0[14] += 1 << 31;
  msg2block1[4] += 1 << 31;
  msg2block1[11] -= 1 << 15;
  msg2block1[14] += 1 << 31;
  if (verbose) cout << endl;
  return true;
}
