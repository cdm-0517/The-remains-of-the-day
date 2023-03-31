## 平面最近点对

```C++
using ll = long long;

ll norm(ll x){
    return x*x;
}

struct point{
    ll x, y;
    bool operator < (const point &t) const{
        return x==t.x?y<t.y:x<t.x;
    }
};

void solve(){
    int n; cin >> n;
    vector<point>vec(n);
    for(auto &[x, y]: vec) cin >> x >> y;
    sort(begin(vec), end(vec));

    ll min_d2 = 1e18, min_d = 1e9;

    queue<point>q;
    set<pair<ll, ll>>now;

    for(auto &[x, y]: vec){
        while(!q.empty() and x - q.front().x > min_d){
            auto [lx, ly] = q.front(); q.pop();
            now.erase({ly, lx});
        }

        auto lit = now.lower_bound({y-min_d-1, 0});
        auto rit = now.lower_bound({y+min_d+1, 0});

        for(auto it = lit ; it != rit ; it++) {
            auto [ty, tx] = *it;

            ll tmp_d2 = norm(x-tx) + norm(y-ty);
            if(tmp_d2 < min_d2){
                min_d2 = tmp_d2;
                min_d = sqrt(tmp_d2 + 1);
            }
        }

        q.push({x, y});
        now.insert({y, x});
    }
    cout << min_d2 << '\n';
}
```

## 最小圆覆盖

```C++
using db = long double;
constexpr db eps = 1e-10;

db norm(db x){
    return x * x;
}

struct point{
    db x, y;
    friend db dis2(const point &A, const point &B){
        return norm(A.x-B.x) + norm(A.y-B.y);
    }
    point operator + (const point &a) const{
        return {x+a.x, y+a.y};
    }
    point operator - (const point &a) const{
        return {x-a.x, y-a.y};
    }
    db operator ^ (const point &a) const{
        return x*a.y - y*a.x;
    }
    point operator * (const db &k) const{
        return {k*x, k*y};
    }
    point rot90() const{
        return {-y, x};
    }
};

struct line {
    point p, v;

    point inter(const line &a){
        return p + v*((a.v^(p-a.p))/(v^a.v));
    }
};

struct circle{
    point o;
    db r;

    circle(const point &O, const db &R):o(O), r(R) {}

    circle (const point &a){
        o = a;
        r = 0;
    }
    circle(const point &a, const point &b){
        o = (a+b) * 0.5;
        r = sqrt(dis2(o, a));
    }
    circle(const point &a, const point &b, const point &c){
        auto A = (a+b) * 0.5;
        auto B = (a+c) * 0.5;

        auto v1 = (a-A).rot90();
        auto v2 = (c-B).rot90();

        o = line{A, v1}.inter(line{B, v2});
        r = sqrt(dis2(o, a));
    }
    circle (vector<point> vec){
        mt19937 rng(114514);
        shuffle(begin(vec), end(vec), rng);

        int n = vec.size();

        *this =circle(vec[0]);
        for(int i = 1 ; i < n ; i ++){
            if((*this).is_in(vec[i])==1) continue;
            *this = circle(vec[i]);
            for(int j = 0 ; j < i ; j ++){
                if((*this).is_in(vec[j])==1) continue;
                *this = circle(vec[i], vec[j]);
                for(int k = 0 ; k < j ; k ++){
                    if((*this).is_in(vec[k])==1) continue;
                    *this = circle(vec[i], vec[j], vec[k]);
                }
            }
        }
    }
    // -1 on, 0 out, 1 in
    int is_in(const point &a) const{
        db d = sqrt(dis2(o, a));
        return abs(d-r)<=eps?-1:d<r-eps;
    }
};
```

## 三维偏序

```C++
// 可能重点，求小于等于的三维偏序
const int maxn = 2e5 + 100;

struct fenwick_tree{
    int val[maxn+100];
    void add(int x, int add){
        for(; x<=maxn; x += (x&-x)) val[x] += add;
    }
    int ask(int x){
        int ans = 0;
        for(;x;x&=x-1) ans += val[x];
        return ans;
    }
    int ask(int l, int r){
        return ask(r) - ask(l-1);
    }
}bit;


void solve(){
    int n, k; cin >> n >> k;
    vector<tuple<int, int, int>>vec(n);
    for(auto &[a, b, c]: vec) cin >> a >> b >> c;

    vector<int>ans(n);

    map<tuple<int, int, int>, int>f;

    function<void(int, int, vector<int>)>func = [&](int l, int r, vector<int>indices){
        if(l==r){
            stable_sort(begin(indices), end(indices), [&](int x, int y){
                const auto &[x1, y1, z1] = vec[x];
                const auto &[x2, y2, z2] = vec[y];
                return tie(y1, z1)<tie(y2, z2);
            });
            for(auto &id: indices){
                const auto &[x, y, z] = vec[id];
                ans[id] += bit.ask(z);
                bit.add(z, 1);
            }
            for(auto &id: indices){
                const auto &[x, y, z] = vec[id];
                bit.add(z, -1);
            }
        }
        else{
            int mid = (l+r)/2;
            vector<int>left_indices, right_indices;
            for(auto &id: indices) {
                const auto &[a, b, c] = vec[id];
                if(a<=mid) left_indices.push_back(id);
                else right_indices.push_back(id);
            }
            func(l, mid, left_indices);
            func(mid+1, r, right_indices);

            stable_sort(begin(left_indices), end(left_indices), [&](int x, int y){
                const auto &[x1, y1, z1] = vec[x];
                const auto &[x2, y2, z2] = vec[y];
                return tie(y1, z1)<tie(y2, z2);
            });
            stable_sort(begin(right_indices), end(right_indices), [&](int x, int y){
                const auto &[x1, y1, z1] = vec[x];
                const auto &[x2, y2, z2] = vec[y];
                return tie(y1, z1)<tie(y2, z2);
            });

            auto tmp = left_indices;
            reverse(begin(tmp), end(tmp));

            vector<int>log;

            for(auto &id: right_indices){
                auto [x1, y1, z1] = vec[id];
                while(tmp.size() and get<1>(vec[tmp.back()])<=y1){
                    const auto [x2, y2, z2] = vec[tmp.back()];
                    bit.add(z2, 1);
                    log.push_back(tmp.back());
                    tmp.pop_back();
                }
                ans[id] += bit.ask(z1);
            }
            for(auto &id: log){
                auto [x, y, z] = vec[id];
                bit.add(z, -1);
            }
        }   
    };
    vector<int>ids(n); iota(begin(ids), end(ids), 0);
    func(1, k, ids);

    for(int i = 0 ; i < n ; i ++){
        f[vec[i]] = max(f[vec[i]], ans[i]);
    }
    vector<int>freq(n);
    for(int i = 0 ; i < n ; i ++){
        freq[f[vec[i]]] ++;
    }
    for(int i = 0 ; i < n ; i ++) cout << freq[i] << '\n';
}
```