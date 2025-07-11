// terminal_histogram.cpp
// A utility for drawing multiple histograms and text boxes in terminal using
// termbox2

// For termbox2: define implementation and compatibility
// #include <__format/format_arg_store.h>
// #define TB_IMPL
// #define TB_OPT_V1_COMPAT
#include <csignal>
#include <string>
#include "ext/termbox2.h"
#include "histogram.h"
// Example
int main() {
    Screen s(2, 2);
    int zo = load_terminfo();
    int h1 = s.addHistogram<int>(0, 0, "Random A"),
        h2 = s.addHistogram<int>(0, 1, "Random B"),
        t = s.addText<int, double, std::string>(1, 0, "Total: {}  Ratio: {}  Note: {}"),
        term = s.addTerm<50>(1, 1, "Terminal");
        
    struct tb_event ev;
    int cnt = 0;
    while (s.running) {
        while (tb_peek_event(&ev, 0) != TB_ERR_NO_EVENT)
            if (ev.type == TB_EVENT_KEY && (ev.key == TB_KEY_CTRL_C || ev.key == TB_KEY_ESC || ev.key == TB_KEY_ENTER) )
                s.running = false;
        s.pushHistogram<int>(h1, rand() % 80);
        s.pushHistogram<int>(h2, rand() % 150);
        s.setTextValues(t, (int)++cnt, (double) 0, std::string{"5"});
        s.printTerminal<50>(term, "Nice !{}\n",cnt);
        s.update();
        ++cnt;
    }
    return 0;
}
