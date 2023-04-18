
#include "States.hpp"

#include <stdexcept>

void lsms::States::relativistic_atomic_states(int core_charge,
                                              std::vector<int> &n,
                                              std::vector<int> &l,
                                              std::vector<int> &spin,
                                              std::vector<int> &kappa,
                                              std::vector<double> &occupation) {
  /*
   * Here the atomic states are calculated in oder to obtain the Dirac
   * Finestructure
   *
   * The degeneracy for n,l quantum level is 2*(2*l+1)
   *
   *
   */

  std::vector<int> n_nonrel;
  std::vector<int> l_nonrel;
  std::vector<double> occupation_nonrel;

  /**
   * First calculate the non relativistic quantum numbers
   */

  n.clear();
  l.clear();
  spin.clear();
  occupation.clear();
  kappa.clear();

  nonrelativistic_atomic_states(core_charge, n_nonrel, l_nonrel,
                                occupation_nonrel);

  double occ_number;

  int number_of_states = n_nonrel.size();

  for (int i = 0; i < number_of_states; ++i) {
    if (l_nonrel[i] == 0) {
      n.emplace_back(n_nonrel[i]);
      l.emplace_back(l_nonrel[i]);
      spin.emplace_back(1);
      occupation.emplace_back(occupation_nonrel[i]);
    } else {
      n.emplace_back(n_nonrel[i]);
      l.emplace_back(l_nonrel[i]);

      /*
       * f_nl = 2*(2*l+1) -> this is due to the magentic quantum numbers
       *
       * In the relativstic formulation there is still a degeneracy
       *
       * f_rel = 2*l/(2*(2*l + 1) + (2*l + 2)/(2*(2*l + 1)
       *
       */

      occ_number = occupation_nonrel[i] * (2 * l_nonrel[i]) /
                   (2 * (2 * l_nonrel[i] + 1));
      occupation.emplace_back(occ_number);
      spin.emplace_back(0);

      n.emplace_back(n_nonrel[i]);
      l.emplace_back(l_nonrel[i]);
      occ_number = occupation_nonrel[i] * (2 * l_nonrel[i] + 2) /
                   (2 * (2 * l_nonrel[i] + 1));
      occupation.emplace_back(occ_number);

      spin.emplace_back(1);
    }
  }

  int number_of_rel_states = n.size();

  for (int i = 0; i < number_of_rel_states; ++i) {
    if (spin[i] == 1) {
      int l_value = -l[i] - 1;
      kappa.emplace_back(l_value);
    } else {
      kappa.emplace_back(l[i]);
    }
  }
}

void lsms::States::nonrelativistic_atomic_states(
    int core_charge, std::vector<int> &n, std::vector<int> &l,
    std::vector<double> &occupation) {
  switch (core_charge) {
    case (1):
      n = {1};
      l = {0};
      occupation = {1};
      break;
    case (2):
      n = {1};
      l = {0};
      occupation = {2};
      break;
    case (3):
      n = {1, 2};
      l = {0, 0};
      occupation = {2, 1};

      break;
    case (4):
      n = {1, 2};
      l = {0, 0};
      occupation = {2, 2};

      break;
    case (5):
      n = {1, 2, 2};
      l = {0, 0, 1};
      occupation = {2, 2, 1};

      break;
    case (6):
      n = {1, 2, 2};
      l = {0, 0, 1};
      occupation = {2, 2, 2};

      break;
    case (7):
      n = {1, 2, 2};
      l = {0, 0, 1};
      occupation = {2, 2, 3};

      break;
    case (8):
      n = {1, 2, 2};
      l = {0, 0, 1};
      occupation = {2, 2, 4};

      break;
    case (9):
      n = {1, 2, 2};
      l = {0, 0, 1};
      occupation = {2, 2, 5};

      break;
    case (10):
      n = {1, 2, 2};
      l = {0, 0, 1};
      occupation = {2, 2, 6};

      break;
    case (11):
      n = {1, 2, 2, 3};
      l = {0, 0, 1, 0};
      occupation = {2, 2, 6, 1};

      break;
    case (12):
      n = {1, 2, 2, 3};
      l = {0, 0, 1, 0};
      occupation = {2, 2, 6, 2};

      break;
    case (13):
      n = {1, 2, 2, 3, 3};
      l = {0, 0, 1, 0, 1};
      occupation = {2, 2, 6, 2, 1};

      break;
    case (14):
      n = {1, 2, 2, 3, 3};
      l = {0, 0, 1, 0, 1};
      occupation = {2, 2, 6, 2, 2};

      break;
    case (15):
      n = {1, 2, 2, 3, 3};
      l = {0, 0, 1, 0, 1};
      occupation = {2, 2, 6, 2, 3};

      break;
    case (16):
      n = {1, 2, 2, 3, 3};
      l = {0, 0, 1, 0, 1};
      occupation = {2, 2, 6, 2, 4};

      break;
    case (17):
      n = {1, 2, 2, 3, 3};
      l = {0, 0, 1, 0, 1};
      occupation = {2, 2, 6, 2, 5};

      break;
    case (18):
      n = {1, 2, 2, 3, 3};
      l = {0, 0, 1, 0, 1};
      occupation = {2, 2, 6, 2, 6};

      break;
    case (19):
      n = {1, 2, 2, 3, 3, 4};
      l = {0, 0, 1, 0, 1, 0};
      occupation = {2, 2, 6, 2, 6, 1};

      break;
    case (20):
      n = {1, 2, 2, 3, 3, 4};
      l = {0, 0, 1, 0, 1, 0};
      occupation = {2, 2, 6, 2, 6, 2};

      break;
    case (21):
      n = {1, 2, 2, 3, 3, 3, 4};
      l = {0, 0, 1, 0, 1, 2, 0};
      occupation = {2, 2, 6, 2, 6, 1, 2};

      break;
    case (22):
      n = {1, 2, 2, 3, 3, 3, 4};
      l = {0, 0, 1, 0, 1, 2, 0};
      occupation = {2, 2, 6, 2, 6, 2, 2};

      break;
    case (23):
      n = {1, 2, 2, 3, 3, 3, 4};
      l = {0, 0, 1, 0, 1, 2, 0};
      occupation = {2, 2, 6, 2, 6, 3, 2};

      break;
    case (24):
      n = {1, 2, 2, 3, 3, 3, 4};
      l = {0, 0, 1, 0, 1, 2, 0};
      occupation = {2, 2, 6, 2, 6, 5, 1};

      break;
    case (25):
      n = {1, 2, 2, 3, 3, 3, 4};
      l = {0, 0, 1, 0, 1, 2, 0};
      occupation = {2, 2, 6, 2, 6, 5, 2};

      break;
    case (26):
      n = {1, 2, 2, 3, 3, 3, 4};
      l = {0, 0, 1, 0, 1, 2, 0};
      occupation = {2, 2, 6, 2, 6, 6, 2};

      break;
    case (27):
      n = {1, 2, 2, 3, 3, 3, 4};
      l = {0, 0, 1, 0, 1, 2, 0};
      occupation = {2, 2, 6, 2, 6, 7, 2};

      break;
    case (28):
      n = {1, 2, 2, 3, 3, 3, 4};
      l = {0, 0, 1, 0, 1, 2, 0};
      occupation = {2, 2, 6, 2, 6, 8, 2};

      break;
    case (29):
      n = {1, 2, 2, 3, 3, 3, 4};
      l = {0, 0, 1, 0, 1, 2, 0};
      occupation = {2, 2, 6, 2, 6, 10, 1};

      break;
    case (30):
      n = {1, 2, 2, 3, 3, 3, 4};
      l = {0, 0, 1, 0, 1, 2, 0};
      occupation = {2, 2, 6, 2, 6, 10, 2};

      break;
    case (31):
      n = {1, 2, 2, 3, 3, 3, 4, 4};
      l = {0, 0, 1, 0, 1, 2, 0, 1};
      occupation = {2, 2, 6, 2, 6, 10, 2, 1};

      break;
    case (32):
      n = {1, 2, 2, 3, 3, 3, 4, 4};
      l = {0, 0, 1, 0, 1, 2, 0, 1};
      occupation = {2, 2, 6, 2, 6, 10, 2, 2};

      break;
    case (33):
      n = {1, 2, 2, 3, 3, 3, 4, 4};
      l = {0, 0, 1, 0, 1, 2, 0, 1};
      occupation = {2, 2, 6, 2, 6, 10, 2, 3};

      break;
    case (34):
      n = {1, 2, 2, 3, 3, 3, 4, 4};
      l = {0, 0, 1, 0, 1, 2, 0, 1};
      occupation = {2, 2, 6, 2, 6, 10, 2, 4};

      break;
    case (35):
      n = {1, 2, 2, 3, 3, 3, 4, 4};
      l = {0, 0, 1, 0, 1, 2, 0, 1};
      occupation = {2, 2, 6, 2, 6, 10, 2, 5};

      break;
    case (36):
      n = {1, 2, 2, 3, 3, 3, 4, 4};
      l = {0, 0, 1, 0, 1, 2, 0, 1};
      occupation = {2, 2, 6, 2, 6, 10, 2, 6};

      break;
    case (37):
      n = {1, 2, 2, 3, 3, 3, 4, 4, 5};
      l = {0, 0, 1, 0, 1, 2, 0, 1, 0};
      occupation = {2, 2, 6, 2, 6, 10, 2, 6, 1};

      break;
    case (38):
      n = {1, 2, 2, 3, 3, 3, 4, 4, 5};
      l = {0, 0, 1, 0, 1, 2, 0, 1, 0};
      occupation = {2, 2, 6, 2, 6, 10, 2, 6, 2};

      break;
    case (39):
      n = {1, 2, 2, 3, 3, 3, 4, 4, 4, 5};
      l = {0, 0, 1, 0, 1, 2, 0, 1, 2, 0};
      occupation = {2, 2, 6, 2, 6, 10, 2, 6, 1, 2};

      break;
    case (40):
      n = {1, 2, 2, 3, 3, 3, 4, 4, 4, 5};
      l = {0, 0, 1, 0, 1, 2, 0, 1, 2, 0};
      occupation = {2, 2, 6, 2, 6, 10, 2, 6, 2, 2};

      break;
    case (41):
      n = {1, 2, 2, 3, 3, 3, 4, 4, 4, 5};
      l = {0, 0, 1, 0, 1, 2, 0, 1, 2, 0};
      occupation = {2, 2, 6, 2, 6, 10, 2, 6, 4, 1};

      break;
    case (42):
      n = {1, 2, 2, 3, 3, 3, 4, 4, 4, 5};
      l = {0, 0, 1, 0, 1, 2, 0, 1, 2, 0};
      occupation = {2, 2, 6, 2, 6, 10, 2, 6, 5, 1};

      break;
    case (43):
      n = {1, 2, 2, 3, 3, 3, 4, 4, 4, 5};
      l = {0, 0, 1, 0, 1, 2, 0, 1, 2, 0};
      occupation = {2, 2, 6, 2, 6, 10, 2, 6, 5, 2};

      break;
    case (44):
      n = {1, 2, 2, 3, 3, 3, 4, 4, 4, 5};
      l = {0, 0, 1, 0, 1, 2, 0, 1, 2, 0};
      occupation = {2, 2, 6, 2, 6, 10, 2, 6, 7, 1};

      break;
    case (45):
      n = {1, 2, 2, 3, 3, 3, 4, 4, 4, 5};
      l = {0, 0, 1, 0, 1, 2, 0, 1, 2, 0};
      occupation = {2, 2, 6, 2, 6, 10, 2, 6, 8, 1};

      break;
    case (46):
      n = {1, 2, 2, 3, 3, 3, 4, 4, 4};
      l = {0, 0, 1, 0, 1, 2, 0, 1, 2};
      occupation = {2, 2, 6, 2, 6, 10, 2, 6, 10};

      break;
    case (47):
      n = {1, 2, 2, 3, 3, 3, 4, 4, 4, 5};
      l = {0, 0, 1, 0, 1, 2, 0, 1, 2, 0};
      occupation = {2, 2, 6, 2, 6, 10, 2, 6, 10, 1};

      break;
    case (48):
      n = {1, 2, 2, 3, 3, 3, 4, 4, 4, 5};
      l = {0, 0, 1, 0, 1, 2, 0, 1, 2, 0};
      occupation = {2, 2, 6, 2, 6, 10, 2, 6, 10, 2};

      break;
    case (49):

      n = {1, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5};
      l = {0, 0, 1, 0, 1, 2, 0, 1, 2, 0, 1};
      occupation = {2, 2, 6, 2, 6, 10, 2, 6, 10, 2, 1};

      break;
    case (50):

      n = {1, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5};
      l = {0, 0, 1, 0, 1, 2, 0, 1, 2, 0, 1};
      occupation = {2, 2, 6, 2, 6, 10, 2, 6, 10, 2, 2};

      break;
    case (51):

      n = {1, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5};
      l = {0, 0, 1, 0, 1, 2, 0, 1, 2, 0, 1};
      occupation = {2, 2, 6, 2, 6, 10, 2, 6, 10, 2, 3};

      break;
    case (52):

      n = {1, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5};
      l = {0, 0, 1, 0, 1, 2, 0, 1, 2, 0, 1};
      occupation = {2, 2, 6, 2, 6, 10, 2, 6, 10, 2, 4};

      break;
    case (53):

      n = {1, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5};
      l = {0, 0, 1, 0, 1, 2, 0, 1, 2, 0, 1};
      occupation = {2, 2, 6, 2, 6, 10, 2, 6, 10, 2, 5};

      break;
    case (54):

      n = {1, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5};
      l = {0, 0, 1, 0, 1, 2, 0, 1, 2, 0, 1};
      occupation = {2, 2, 6, 2, 6, 10, 2, 6, 10, 2, 6};

      break;
    case (55):

      n = {1, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 6};
      l = {0, 0, 1, 0, 1, 2, 0, 1, 2, 0, 1, 0};
      occupation = {2, 2, 6, 2, 6, 10, 2, 6, 10, 2, 6, 1};

      break;
    case (56):

      n = {1, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 6};
      l = {0, 0, 1, 0, 1, 2, 0, 1, 2, 0, 1, 0};
      occupation = {2, 2, 6, 2, 6, 10, 2, 6, 10, 2, 6, 2};

      break;
    case (57):

      n = {1, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5, 6};
      l = {0, 0, 1, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0};
      occupation = {2, 2, 6, 2, 6, 10, 2, 6, 10, 2, 6, 1, 2};

      break;
    case (58):

      n = {1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 6};
      l = {0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 0};
      occupation = {2, 2, 6, 2, 6, 10, 2, 6, 10, 1, 2, 6, 1, 2};

      break;
    case (59):

      n = {1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 6};
      l = {0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 0};
      occupation = {2, 2, 6, 2, 6, 10, 2, 6, 10, 3, 2, 6, 2};

      break;
    case (60):

      n = {1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 6};
      l = {0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 0};
      occupation = {2, 2, 6, 2, 6, 10, 2, 6, 10, 4, 2, 6, 2};

      break;
    case (61):

      n = {1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 6};
      l = {0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 0};
      occupation = {2, 2, 6, 2, 6, 10, 2, 6, 10, 5, 2, 6, 2};

      break;
    case (62):

      n = {1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 6};
      l = {0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 0};
      occupation = {2, 2, 6, 2, 6, 10, 2, 6, 10, 6, 2, 6, 2};

      break;
    case (63):

      n = {1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 6};
      l = {0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 0};
      occupation = {2, 2, 6, 2, 6, 10, 2, 6, 10, 7, 2, 6, 2};

      break;
    case (64):

      n = {1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 6};
      l = {0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 0};
      occupation = {2, 2, 6, 2, 6, 10, 2, 6, 10, 7, 2, 6, 1, 2};

      break;
    case (65):

      n = {1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 6};
      l = {0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 0};
      occupation = {2, 2, 6, 2, 6, 10, 2, 6, 10, 9, 2, 6, 2};

      break;
    case (66):

      n = {1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 6};
      l = {0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 0};
      occupation = {2, 2, 6, 2, 6, 10, 2, 6, 10, 10, 2, 6, 2};

      break;
    case (67):

      n = {1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 6};
      l = {0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 0};
      occupation = {2, 2, 6, 2, 6, 10, 2, 6, 10, 11, 2, 6, 2};

      break;
    case (68):

      n = {1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 6};
      l = {0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 0};
      occupation = {2, 2, 6, 2, 6, 10, 2, 6, 10, 12, 2, 6, 2};

      break;
    case (69):

      n = {1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 6};
      l = {0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 0};
      occupation = {2, 2, 6, 2, 6, 10, 2, 6, 10, 13, 2, 6, 2};

      break;
    case (70):

      n = {1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 6};
      l = {0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 0};
      occupation = {2, 2, 6, 2, 6, 10, 2, 6, 10, 14, 2, 6, 2};

      break;
    case (71):

      n = {1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 6};
      l = {0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 0};
      occupation = {2, 2, 6, 2, 6, 10, 2, 6, 10, 14, 2, 6, 1, 2};

      break;
    case (72):

      n = {1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 6};
      l = {0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 0};
      occupation = {2, 2, 6, 2, 6, 10, 2, 6, 10, 14, 2, 6, 2, 2};

      break;
    case (73):

      n = {1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 6};
      l = {0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 0};
      occupation = {2, 2, 6, 2, 6, 10, 2, 6, 10, 14, 2, 6, 3, 2};

      break;
    case (74):

      n = {1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 6};
      l = {0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 0};
      occupation = {2, 2, 6, 2, 6, 10, 2, 6, 10, 14, 2, 6, 4, 2};

      break;
    case (75):

      n = {1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 6};
      l = {0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 0};
      occupation = {2, 2, 6, 2, 6, 10, 2, 6, 10, 14, 2, 6, 5, 2};

      break;
    case (76):

      n = {1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 6};
      l = {0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 0};
      occupation = {2, 2, 6, 2, 6, 10, 2, 6, 10, 14, 2, 6, 6, 2};

      break;
    case (77):

      n = {1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 6};
      l = {0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 0};
      occupation = {2, 2, 6, 2, 6, 10, 2, 6, 10, 14, 2, 6, 7, 2};

      break;
    case (78):

      n = {1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 6};
      l = {0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 0};
      occupation = {2, 2, 6, 2, 6, 10, 2, 6, 10, 14, 2, 6, 9, 1};

      break;
    case (79):

      n = {1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 6};
      l = {0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 0};
      occupation = {2, 2, 6, 2, 6, 10, 2, 6, 10, 14, 2, 6, 10, 1};

      break;
    case (80):

      n = {1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 6};
      l = {0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 0};
      occupation = {2, 2, 6, 2, 6, 10, 2, 6, 10, 14, 2, 6, 10, 2};

      break;
    case (81):

      n = {1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 6, 6};
      l = {0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 0, 1};
      occupation = {2, 2, 6, 2, 6, 10, 2, 6, 10, 14, 2, 6, 10, 2, 1};

      break;
    case (82):

      n = {1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 6, 6};
      l = {0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 0, 1};
      occupation = {2, 2, 6, 2, 6, 10, 2, 6, 10, 14, 2, 6, 10, 2, 2};

      break;
    case (83):

      n = {1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 6, 6};
      l = {0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 0, 1};
      occupation = {2, 2, 6, 2, 6, 10, 2, 6, 10, 14, 2, 6, 10, 2, 3};

      break;
    case (84):

      n = {1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 6, 6};
      l = {0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 0, 1};
      occupation = {2, 2, 6, 2, 6, 10, 2, 6, 10, 14, 2, 6, 10, 2, 4};

      break;
    case (85):

      n = {1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 6, 6};
      l = {0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 0, 1};
      occupation = {2, 2, 6, 2, 6, 10, 2, 6, 10, 14, 2, 6, 10, 2, 5};

      break;
    case (86):

      n = {1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 6, 6};
      l = {0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 0, 1};
      occupation = {2, 2, 6, 2, 6, 10, 2, 6, 10, 14, 2, 6, 10, 2, 6};

      break;
    case (87):

      n = {1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 6, 6, 7};
      l = {0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 0, 1, 0};
      occupation = {2, 2, 6, 2, 6, 10, 2, 6, 10, 14, 2, 6, 10, 2, 6, 1};

      break;
    case (88):

      n = {1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 6, 6, 7};
      l = {0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 0, 1, 0};
      occupation = {2, 2, 6, 2, 6, 10, 2, 6, 10, 14, 2, 6, 10, 2, 6, 2};

      break;
    case (89):

      n = {1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 6, 6, 6, 7};
      l = {0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 0, 1, 2, 0};
      occupation = {2, 2, 6, 2, 6, 10, 2, 6, 10, 14, 2, 6, 10, 2, 6, 1, 2};

      break;
    case (90):

      n = {1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 6, 6, 6, 7};
      l = {0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 0, 1, 2, 0};
      occupation = {2, 2, 6, 2, 6, 10, 2, 6, 10, 14, 2, 6, 10, 2, 6, 2, 2};

      break;
    case (91):

      n = {1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 6, 6, 6, 7};
      l = {0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 0};
      occupation = {2, 2, 6, 2, 6, 10, 2, 6, 10, 14, 2, 6, 10, 2, 2, 6, 1, 2};

      break;
    case (92):

      n = {1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 6, 6, 6, 7};
      l = {0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 0};
      occupation = {2, 2, 6, 2, 6, 10, 2, 6, 10, 14, 2, 6, 10, 3, 2, 6, 1, 2};
      break;

    default:
      throw std::invalid_argument("This core charge number is not supported");
  }
}
