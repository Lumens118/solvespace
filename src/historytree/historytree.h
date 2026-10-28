#include "sketch.h"

class HistMgr {
public:
    static void MoveUp(hGroup hg);
    static void MoveDown(hGroup hg);

private:
    static void Move(hGroup hg, MoveDir dir);  // 唯一核心
    static const char *MoveErrText(MoveErr e); // 唯一的错误码→文案映射
};