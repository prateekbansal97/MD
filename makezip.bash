rm -r zip
mkdir zip
cd zip
for i in `ls ../*.cpp` ; do cp $i . ; done
for i in `ls ../*.h` ; do cp $i . ; done
cd ..
tar -czf ./md.tar.gz ./zip
