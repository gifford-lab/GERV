// Copyright 2014 Daniel Kang

#include <string.h>
#include <stdio.h>
#include <ctype.h>

char table[256];
char table2[256];

char rep(int c) {
  c = toupper(c);
  c = table[c];
  c = table2[c];
  return c;
}

// I apologize to the gods of programming for the abuse of global variables
void set_table() {
  int i;
  for (i = 0; i < 256; i++)
    table[i] = 'N';
  table['A'] = 'A';
  table['T'] = 'T';
  table['C'] = 'C';
  table['G'] = 'G';

  for (i = 0; i < 256; i++)
    table2[i] = 4;
  table2['A'] = 0;
  table2['T'] = 1;
  table2['G'] = 2;
  table2['C'] = 3;
}

void proc_chr(FILE *fout, char *fin_name) {
  int size;
  char buf[2000];

  FILE *fin = fopen(fin_name, "r");

  while (!feof(fin)) {
    memset(buf, 0, sizeof(buf));
    if (fgets(buf, sizeof(buf), fin) == NULL) break;
    if (buf[0] == '>') continue;

    size = strlen(buf);
    if (buf[size-1] == '\n') size--;

    for (int i = 0; i < size; i++)
      buf[i] = rep(buf[i]);

    fwrite(buf, 1, size, fout);
  }

  fclose(fin);
}

int main(int argc, char **argv) {
  FILE *fin, *fout;
  char buf[2000];
  int size;

  if (argc != 3) {
    fprintf(stdout, "Usage: %s <out> <in>\n", argv[0]);
    return 0;
  } else {
    fout = fopen(argv[1], "w");
    fin = fopen(argv[2], "r");
  }

  set_table();

  while (!feof(fin)) {
    memset(buf, 0, 2000);
    if (fgets(buf, 1024, fin) == NULL) break;

    size = strlen(buf);
    if (buf[size-1] == '\n') buf[size-1] = 0;
    if (size < 2) break;

    fprintf(stdout, "Processing: %s\n", buf);
    proc_chr(fout, buf);

    if (feof(fin)) break;
  }

  fclose(fin);
  fclose(fout);

  return 0;
}

