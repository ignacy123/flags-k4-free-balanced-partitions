from flag_utils import read_flags
from argparse import ArgumentParser

f_expression = "result.txt"
f_values = "tripartite.txt"


def main():
    parser = ArgumentParser()
    parser.add_argument("--values", type=str, required=True)
    parser.add_argument("--expression", type=str, required=True)

    args = parser.parse_args()

    with open(args.expression) as f1:
        with open(args.values) as f2:
            sum = 0
            for flag1, flag2 in zip(read_flags(f1), read_flags(f2)):
                if not flag1.edges == flag2.edges:
                    print("AAA")
                val = flag1.prob * flag2.prob
                sum += val
                print(flag1.prob, flag2.prob)
                # if abs(val) > 0.0001:
                #     print(flag1)
                #     print(val)

        print(sum)


if __name__ == "__main__":
    main()
