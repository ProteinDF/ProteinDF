#include "PdfUtils.h"
#include "TlUtils.h"

bool PdfUtils::isComment(const std::string& str)
{
    std::string check = str;
    const int length = str.size();

    TlUtils::trim_ws(check);
    if ((length >= 1) && (check.substr(0, 1) == "#")) {
        return true;
    }
    if (length >= 2) {
        const std::string head = check.substr(0, 2);
        if ((head == "//") || (head == "--")) {
            return true;
        }
    }

    return false;
}
