#include "DfInputdata_Parallel.h"
#include "TlCommunicate.h"

DfInputdata_Parallel::DfInputdata_Parallel()
{
}


DfInputdata_Parallel::~DfInputdata_Parallel()
{
}


TlSerializeData DfInputdata_Parallel::main()
{
    TlCommunicate& rComm = TlCommunicate::getInstance();

    if (rComm.isMaster() == true) {
        // masterのthis->data_にデータが格納される
        DfInputdata::main();
    }

    rComm.broadcast(this->data_);

    return this->data_;
}


