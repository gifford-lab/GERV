import argparse, csv
import numpy as np

def parse_command_line():
  parser = argparse.ArgumentParser(description='')
  parser.add_argument('--raw_file', required=True, help='Raw input')
  parser.add_argument('--lengths', required=True, help='lengths of strings')
  parser.add_argument('--csv_file', required=True, help='CSV output')
  parser.add_argument('--padding', required=True, help='Padding to ignore')
  args = parser.parse_args()
  return args

def main():
  args = parse_command_line()

  padding = int(args.padding)
  arr = np.fromfile(args.raw_file, dtype='float32')
  arr = arr[padding:-padding]

  writer = csv.writer(open(args.csv_file, 'wb'), delimiter=',')
  writer.writerow(['String number', 'Position', 'Score'])


  offset = 0
  for str_num, line in enumerate(open(args.lengths)):
    str_len = int(line.strip())

    for i in range(str_len):
      writer.writerow([str_num, i, arr[offset + i]])

    offset += padding + str_len


if __name__ == '__main__':
  main()


