    static inline int64 good_size(const int64 need) {
        int64 good = 1;
        while (good < need) good *= 2;
        if (good < 8) return good;
        if (good*5/8 >= need) return good*5/8;
        if (good*3/4 >= need) return good*3/4;
        return good;
    }
