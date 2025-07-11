// terminal_histogram.cpp
// A utility for drawing multiple histograms and text boxes in terminal using
// termbox2

// For termbox2: define implementation and compatibility
// #include <__format/format_functions.h>
#include <format>
// #include <print>
#define TB_IMPL
#define TB_OPT_V1_COMPAT
#include <algorithm>
#include <csignal>
#include <cstdint>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include "circularqueue.h"
#include "ext/termbox2.h"
class Cell {
   public:
    int gr, gc;  // grid coords
    Cell(int r, int c) : gr(r), gc(c) {}
    virtual ~Cell() = default;
    virtual void draw(int top, int left, int height, int width) = 0;
};

template <typename T>
const T &unmove(T &&x) {
    return x;
}

// TextCell: uses std::format for placeholders `{}` and tb_printf for rendering
template <typename... Args>
class TextCell : public Cell {
    std::string fmt;

    std::tuple<Args...> store;

   public:
    TextCell(int r, int c, const std::string &formatStr)
        : Cell(r, c), fmt(formatStr), store(Args{}...) {}

    void setValues(const Args &...args) {
        // vals.clear();
        store = std::make_tuple(args...);
        // (vals.emplace_back(std::forward<Args>(args)), ...);
    }

    static auto my_mfa(Args const &...args) {
        return std::make_format_args(args...);
    }
    void draw(int top, int left, int height, int width) override {
        // prepare arg store

        std::string text = std::vformat(fmt, std::apply(my_mfa, store));
        // wrap and print

        std::stringstream ss{text};
        int y = 0;
        std::string line;
        while (getline(ss, line, '\n')) {
            int pos = y * width;
            // std::string line = (pos < (int)text.size())
            //                        ? text.substr(pos, width)
            //                        : std::string();
            tb_printf(left + 2, top + y + 2, TB_WHITE, TB_DEFAULT, "%s",
                      line.c_str());
            ++y;
        }
    }
};

template<int HL>
class TerminalCell : public Cell {
    std::string name;
    cdqa_t<std::string, HL> que;
   public:
    TerminalCell(int r, int c, const std::string &name)
        : Cell(r, c), name{name}{}

    template<class... Args>
    void print(std::format_string<Args...> fmt, Args &&... args){
        que.push_front(std::format(fmt, args...));
        if(que.size() == HL){
            que.pop_back();
        }
    }
    void draw(int top, int left, int height, int width) override {
        int i = 0;
        int idx = 0;
        while(i < height){
            for(int j = 0; j < que[idx].size(); j+=width){
                tb_print(left+1, top+height-i-j+1, TB_WHITE, TB_DEFAULT, que[idx].substr(j,width).c_str());
                ++i;
                if(i >= height){ break;}
            }
            ++idx;   
        }
    }
};

// HistogramCell: integer histogram with title, X/Y axes
template <typename T>
class HistogramCell : public Cell {
    std::string title;
    std::vector<T> data;

   public:
    HistogramCell(int r, int c, const std::string &t) : Cell(r, c), title(t) {}
    void push(T v) { data.push_back(v); }

    void draw(int top, int left, int height, int width) override {
        if (data.empty()) return;
        // compute range and bins
        auto [minIt, maxIt] = std::minmax_element(data.begin(), data.end());
        T minv = *minIt, maxv = *maxIt;
        if (minv == maxv) maxv = minv + 1;
        int bins = width - 5;
        if (bins < 1) return;
        std::vector<int> cnt(bins);
        for (auto &v : data) {
            int idx = int((double)(v - minv) / (maxv - minv) * (bins - 1));
            cnt[idx]++;
        }
        int maxc = *std::max_element(cnt.begin(), cnt.end());
        if (maxc == 0) maxc = 1;
        int barH = height - 2;
        int br = top + 1, bc = left + 5;
        // bars
        for (int x = 0; x < bins; ++x) {
            int h = int((double)cnt[x] / maxc * barH);
            for (int y = 0; y < barH; ++y) {
                uint32_t ch = (y >= barH - h) ? U'█' : ' ';
                tb_set_cell(bc + x, br + y, ch, TB_YELLOW, TB_DEFAULT);
            }
        }
        // X labels
        int stepX = std::max(1, bins / 5);
        for (int x = 0; x < bins; x += stepX) {
            T v = T(minv + (double)(x + 1) / bins * (maxv - minv));
            tb_printf(bc + x + 1, br + barH + 1, TB_CYAN, TB_DEFAULT, "%d",
                      (int)v);
        }
        // Y axis
        int steps = 5;
        for (int i = 0; i <= steps; ++i) {
            int y = br + barH - (barH * i / steps);
            int v = int(maxc * i / steps);
            tb_printf(left + 1, y + 1, TB_CYAN, TB_DEFAULT, "%d", v);
            tb_set_cell(left + 4, y, U'│', TB_WHITE, TB_DEFAULT);
        }  // title
        int tx = left + (width - (int)title.size()) / 2;
        tb_printf(tx + 1, top + 1, TB_RED, TB_DEFAULT, "%s", title.c_str());
    }
};

// Screen manager
class Screen {
    int gridR, gridC;
    int cH, cW;
    std::vector<std::shared_ptr<Cell>> cells;

   public:
    bool running = true;
    Screen(int r, int c) : gridR(r), gridC(c) {
        if (tb_init()) throw;
        
        tb_clear();
        std::signal(SIGINT, [](int) {
            tb_shutdown();
            std::exit(0);
        });
        tb_set_input_mode(TB_INPUT_MOUSE );

        compute();
    }
    ~Screen() { tb_shutdown(); }
    template <typename T>
    int addHistogram(int gr, int gc, const std::string &t) {
        cells.push_back(std::make_shared<HistogramCell<T>>(gr, gc, t));
        return cells.size() - 1;
    }
    template <class... Args>
    int addText(int gr, int gc, const std::string &f) {
        cells.push_back(std::make_shared<TextCell<Args...>>(gr, gc, f));
        return cells.size() - 1;
    }

    template<int N>
    int addTerm(int gr, int gc, const std::string &n){
        cells.push_back(std::make_shared<TerminalCell<N>>(gr,gc,n));
        return cells.size() - 1;
    }

    template <typename... Args>
    void setTextValues(int id, Args... args) {
        auto t = dynamic_cast<TextCell<Args...> *>(cells[id].get());
        if (!t) throw;
        t->setValues(std::forward<Args>(args)...);
    }
    template <typename T>
    void pushHistogram(int id, T v) {
        auto h = dynamic_cast<HistogramCell<T> *>(cells[id].get());
        if (!h) throw;
        h->push(v);
    }
    template <int N, class ...Args>
    void printTerminal(int id, std::format_string<Args... > fmt, Args &&...args) {
        auto h = dynamic_cast<TerminalCell<N> *>(cells[id].get());
        if (!h) throw;
        h->print(fmt,args...);
    }
    

    void update() {
        compute();
        tb_clear();
        for (auto &c : cells) {
            int top = c->gr * cH, left = c->gc * cW;
            int h = cH - (c->gr < gridR - 1), w = cW - (c->gc < gridC - 1);
            c->draw(top, left, h, w);
        }
        drawBorders();
        tb_present();
    }

   private:
    void compute() {
        cW = tb_width() / gridC;
        cH = tb_height() / gridR;
    }
    void drawBorders() {
        int TW = tb_width(), TH = tb_height();
        for (int r = 1; r < gridR; ++r) {
            int y = r * cH;
            for (int x = 0; x < TW; ++x)
                tb_set_cell(x, y, U'─', TB_WHITE, TB_DEFAULT);
        }
        for (int c = 1; c < gridC; ++c) {
            int x = c * cW;
            for (int y = 0; y < TH; ++y)
                tb_set_cell(x, y, U'│', TB_WHITE, TB_DEFAULT);
        }
        for (int r = 1; r < gridR; ++r)
            for (int c = 1; c < gridC; ++c)
                tb_set_cell(c * cW, r * cH, U'┼', TB_WHITE, TB_DEFAULT);
    }
};