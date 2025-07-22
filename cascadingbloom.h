#ifndef _CASCADING_BLOOM_H
#define _CASCADING_BLOOM_H


template<class first, class ...BF>
struct cascadebloom :  cascadebloom<BF...>{
    first bf;
    cascadebloom<first, BF...> (first f, BF... b) : cascadebloom<BF...>(b...), bf(f){}

    bool contains(const auto &element){
        if (bf.contains(element)){
            return cascadebloom<BF...>::contains(element);
        }
        return false;
    }

    bool insert(const auto &element){
        if(bf.contains(element)){
            return cascadebloom<BF...>::insert(element);
        }
        bf.insert(element);
        return false;
    }
};

template<class BF>
struct cascadebloom<BF>{
    BF bf;
    cascadebloom(BF b) : bf(b) {}
    bool contains(const auto &element){
        return bf.contains(element);
    }
    bool insert(const auto &element){
        if(bf.contains(element)){
            return true;
        }
        bf.insert(element);
        return false;
    }
};
#endif