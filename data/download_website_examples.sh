
for name in bv_10k\
            bv_100k\
            bv_1M\
            tree_d9_c4\
            tree_d6_c5-10\
            tree_d8_c5;
do
  curl --request GET -sL \
       --url "https://algo2.iti.kit.edu/download/example_${name}.txt"\
       --output "example_${name}.txt"
done