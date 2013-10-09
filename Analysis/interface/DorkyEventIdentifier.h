#ifndef DORKYEVENTIDENTFIER_H
#define DORKYEVENTIDENTFIER_H

// this is a workaround for removing duplicate events 

struct DorkyEventIdentifier 
{
    unsigned long int run, event, lumi;
    bool operator < (const DorkyEventIdentifier &) const;
    bool operator == (const DorkyEventIdentifier &) const;
};

bool is_duplicate(const DorkyEventIdentifier &id);

void reset();

#endif // AT_DORKYEVENTIDENTFIER_H
