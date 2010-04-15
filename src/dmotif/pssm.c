/* $Id: pssm.c,v 1.2 2009/02/19 19:38:48 agd27 Exp $ */

#include <pssm.h>
#include <msa.h>
#include <stringsplus.h>
#include <lists.h>
#include <prob_vector.h>

/* if alphabet == NULL, default will be used */
PSSM *mot_new(int width, char *alphabet) {
  int i;
  PSSM *m = (PSSM*)smalloc(sizeof(PSSM));
  if (alphabet == NULL)
    die("Alphabet cannot be null in mot_new.\n");
  m->width = width;
  m->alphabet = (char*)smalloc((strlen(alphabet) + 1)* sizeof(char));
  m->alphabet[strlen(alphabet)] = '\0';
  strncpy(m->alphabet, alphabet, strlen(alphabet));
  m->alphsize = strlen(m->alphabet);
  m->probs = smalloc(m->width * sizeof(void*));
  for (i = 0; i < m->width; i++)
    m->probs[i] = vec_new(m->alphsize);
  return m;
}

PSSM *mot_read(FILE *F) {
  Regex *width_re = str_re_new("#[[:space:]]*WIDTH[[:space:]]*=[[:space:]]*([0-9]+)");
  Regex *alph_re = str_re_new("#[[:space:]]*ALPHABET[[:space:]]*=[[:space:]]*([A-Z]+)");
  String *line = str_new(STR_MED_LEN), *alphabet = NULL;
  int width = -1, pos = 0, i;
  List *matches = lst_new_ptr(4);
  PSSM *m = NULL;

  while (str_readline(line, F) != EOF) {
    str_trim(line);
    if (line->length == 0) continue;

    if (m == NULL) {
      if (str_re_match(line, width_re, matches, 1) >= 0) 
        str_as_int(lst_get_ptr(matches, 1), &width);
      else if (str_re_match(line, alph_re, matches, 1) >= 0) 
        alphabet = str_dup(lst_get_ptr(matches, 1));
      else 
        die("ERROR: bad header in motif file.\n");

      if (width > 0 && alphabet != NULL)
        m = mot_new(width, alphabet->chars);
    }
    else {
      if (pos >= width)
        die("ERROR: motif probabilities exceed motif width.\n");

      if (str_split(line, NULL, matches) != m->alphsize)
        die("ERROR: number of probabilities per line does not equal alphabet size.\n");

      for (i = 0; i < m->alphsize; i++) {
        if (str_as_dbl(lst_get_ptr(matches, i), &m->probs[pos]->data[i]) != 0)
          die("ERROR: bad symbol ('%s')\n", 
              ((String*)lst_get_ptr(matches, i))->chars);
      }

      pv_normalize(m->probs[pos]);
      pos++;
    }
    lst_free_strings(matches);
    lst_clear(matches);      
  }
  if (m == NULL)
    die("ERROR: incomplete header in motif file.\n");
  if (pos < width)
    die("ERROR: too few motif probability lines.\n");

  str_re_free(width_re);
  str_re_free(alph_re);
  str_free(line);
  str_free(alphabet);
  lst_free_strings(matches);
  lst_free(matches);
  return(m);
}

void mot_write(FILE *F, PSSM *m) {
  int i, j;
  fprintf(F, "# WIDTH = %d\n", m->width);
  fprintf(F, "# ALPHABET = %s\n", m->alphabet);
  for (i = 0; i < m->width; i++) {
    for (j = 0; j < m->alphsize; j++) 
      fprintf(F, "%f%c", m->probs[i]->data[j], j < m->alphsize - 1 ? '\t' : '\n');
  }
}

void mot_free(PSSM *m) {
  int i;
  for (i = 0; i < m->width; i++)
    vec_free(m->probs[i]);
  free(m->probs);
  free(m->alphabet);
  free(m);
}
