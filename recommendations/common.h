 #define DIM 4

 static constexpr size_t precision = 16;// 64 / 3 - std::log2(DIM) - 0.5;
 const size_t nusers   = 943;
 const size_t nitems   = 1ULL << 10;
 const size_t nqueries = 1ULL << 12;

constexpr size_t max_norm = 1ULL << (63 - (3 * precision));
constexpr size_t bucket_size_ = std::ceil(double(max_norm/1.1));
constexpr size_t n_buckets = std::ceil(double(nqueries)/bucket_size_);
constexpr size_t last_buck_size = 0;//nqueries - (n_buckets - 1) * bucket_size_;
constexpr size_t total_bucket_elements =  n_buckets * nitems ; //  bucket_size_ * (n_buckets - 1) + last_buck_size;
constexpr size_t n_cores = 1; 
constexpr size_t thread_cap = std::ceil(double(nqueries)/n_cores); 
constexpr size_t n_threads = nqueries/thread_cap;
constexpr size_t l_thread_cap = nqueries - (thread_cap * (n_threads - 1)); 

  
