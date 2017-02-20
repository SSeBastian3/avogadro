/******************************************************************************

  This source file is part of the Avogadro project.

  Copyright 2016 Kitware, Inc.

  This source code is released under the New BSD License, (the "License").

  Unless required by applicable law or agreed to in writing, software
  distributed under the License is distributed on an "AS IS" BASIS,
  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
  See the License for the specific language governing permissions and
  limitations under the License.

******************************************************************************/

#ifndef GUESSPROJECTIONS_H
#define GUESSPROJECTIONS_H

inline static bool alreadyUsed(const std::string& symbol,
                               const std::vector<std::string>& alreadyLookedAt)
{
  for (size_t i = 0; i < alreadyLookedAt.size(); ++i) {
    if (symbol == alreadyLookedAt[i])
      return true;
  }
  return false;
}

/* Returns a guess for the yaehmop projections QString of a set of atomic
 * symbols. Does them by types. For example, {O, Ti, O} would return:
 * # O
 * atom 1 1.0, 3 1.0
 * # Ti
 * atom 2 1.0
 *
 * @param atomicSymbols The vector of atomic symbols.
 * @return The guessed atom projections. We assume summing (i. e., all
 *         weights are 1.0).
 */
static QString guessTypedAtomProjections(
                           const std::vector<std::string>& atomicSymbols)
{
  QString ret;
  std::vector<std::string> alreadyLookedAt;
  for (size_t i = 0; i < atomicSymbols.size(); ++i) {
    if (!alreadyUsed(atomicSymbols[i], alreadyLookedAt)) {
      ret += (QString("# ") + atomicSymbols[i].c_str() + "\n");
      ret += (QString("atom ") + QString::number(i + 1) + " 1.0");
      for (size_t j = i + 1; j < atomicSymbols.size(); ++j) {
        if (atomicSymbols[i] == atomicSymbols[j])
          ret += (QString(", ") + QString::number(j + 1) + " 1.0");
      }
      ret += "\n";
      alreadyLookedAt.push_back(atomicSymbols[i]);
    }
  }
  return ret;
}

#endif
