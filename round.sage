from argparse import ArgumentParser

def main():
    parser = ArgumentParser()
    parser.add_argument("-s", "--sdp-directory", type=str, required=True)
    parser.add_argument("-f", "--scaling-factor", type=int, required=True)

    args = parser.parse_args()

    load("rounding_Integer.sage")
    round_program(f"SDP/{args.sdp_directory}/problem.dat-s", args.scaling_factor, False)


if __name__ == "__main__":
    main()
