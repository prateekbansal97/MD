rm -r zip
mkdir zip
cd zip
for i in `ls ../src/AmberTopology/*.cpp` ; do cp $i . ; done
for i in `ls ../include/AmberTopology/*.h` ; do cp $i . ; done
cd ..
tar -czf ./md.tar.gz ./zip
