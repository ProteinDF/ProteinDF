#!/bin/sh

dim_list="1000 2000 4000 6000 8000 10000 20000"

echo ">>>> normal matrix bench (on memory)"
for dim in $dim_list; do
    /usr/bin/time ./mmap_bench -d $dim
    echo
done
echo


echo ">>>> symmetric matrix bench (on memory)"
for dim in $dim_list; do
    /usr/bin/time ./mmap_bench -d $dim -s
    echo
done
echo


echo ">>>> normal matrix bench (on disk)"
for dim in $dim_list; do
    /usr/bin/time ./mmap_bench -d $dim -m
    rm /tmp/mmap*
    echo
done
echo


echo ">>>> symmetric matrix bench (on disk)"
for dim in $dim_list; do
    /usr/bin/time ./mmap_bench -d $dim -s -m
    rm /tmp/mmap*
    echo
done
echo


echo ">>>> normal matrix bench (on memory, ramdom) "
for dim in $dim_list; do
    /usr/bin/time ./mmap_bench -d $dim
    echo
done
echo


echo ">>>> symmetric matrix bench (on memory, random)"
for dim in $dim_list; do
    /usr/bin/time ./mmap_bench -d $dim -s
    echo
done
echo


echo ">>>> normal matrix bench (on disk, random)"
for dim in $dim_list; do
    /usr/bin/time ./mmap_bench -d $dim -m
    rm /tmp/mmap*
    echo
done
echo


echo ">>>> symmetric matrix bench (on disk, random)"
for dim in $dim_list; do
    /usr/bin/time ./mmap_bench -d $dim -s -m
    rm /tmp/mmap*
    echo
done
echo

