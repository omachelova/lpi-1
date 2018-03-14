import os.path
import sys
sys.path[0:0] = [os.path.join(sys.path[0], '../../examples/sat')]

import sat


class SudokuSolver(object):
    def __init__(self):
        self.N = 0

    def q(self, n,r,c):
        #vsetko budu uniformne 3ciferne cisla.
        #pracujeme s polom 9*9, lahsie sa to bude citat
        return r *100 + c*10 + n;

    def solve(self, pole):
        self.pole = pole
        solver = sat.SatSolver()
        w = sat.DimacsWriter('sudokusolver_cnf_in.txt')

        for i in range(9):
            for j in range (9):
                if (self.pole[i][j]>0):
                    w.writeLiteral(self.q(int(self.pole[i][j]),i+1,j+1))
                    w.finishClause()
        
        
        # pre cisla 1 az 9
        for n in range (1,10): 
            # v kazdom riadku
            for r in range(1,10):
                # je aspon jedno cislo n
                for c in range(1,10):
                    w.writeLiteral(self.q(n,r,c))
                w.finishClause()

        # pre cisla 1 az 9
        for n in range (1,10): 
            # v kazdom riadku
            for r in range(1,10):
                # na dvoch roznych poziciach
                for c1 in range(1,10):
                    for c2 in range(c1):
                        # nie je na oboch rovnake cislo
                        # q(n,r,c1) => -q(n,r,c2)
                        w.writeImpl(self.q(n,r,c1), -self.q(n,r,c2))
                        
        # pre cisla 1 az 9
        for n in range (1,10): 
            # v kazdom stlpci
            for c in range(1,10):
                # na dvoch roznych miestach
                for r1 in range(1,10):
                    for r2 in range(r1):
                        # nie je na oboch rovnake cislo
                        w.writeImpl(self.q(n,r1,c), -self.q(n,r2,c))

        # uhlopriecky
        # trochu neefektivnejsie, ale rychlejsie na napisanie
#        for c1 in range(N):
#            for c2 in range(N):
#                for r1 in range(N):
#                    for r2 in range(N):
                        # rozne pozicie
#                        if (self.q(r1,c1) != self.q(r2,c2)):
                            # r1,c1 a r2,c2 su na uhlopriecke
#                            if (r1+c1 == r2+c2) or (r1+c2 == r2+c1):
#                                w.writeImpl(self.q(r1,c1), -self.q(r2,c2))

        #podstvorce
        #pre kazde cislo
        for n in range (1,10): 
            #su 3 rady a 3 stlpce podstvorcov
            for i in range (3):
                for j in range (3):
                    #kazdy podstvorec je 3x3
                    for k in range (3):
                        for l in range (3):
                            #stlpec/riadok sa da vypocitat ako podstvorec*3+pozicia v podstvorci+1
                            w.writeLiteral(self.q(n,r,c))
                    w.finishClause()



        w.close()
        ok, sol = solver.solve(w, 'sudoku_out.txt')

        ret = [[]]
        if ok:
            for x in sol:
                if x>0:
                    x -= 1
                    ret.append( (x // 10, x % 10, x // 100, x % 100) )
        return ret


if __name__ == '__main__':
    Pole=[[0,0,0,0,0,0,0,0,0],
          [0,0,0,0,0,0,0,0,0],
          [0,0,0,0,0,0,0,0,0],
          [0,0,0,0,0,0,0,0,0],
          [0,0,0,0,0,0,0,0,0],
          [0,0,0,0,0,0,0,0,0],
          [0,0,0,0,0,0,0,0,0],
          [0,0,0,0,0,0,0,0,0],
          [0,0,0,0,0,0,0,0,0]]
    
    sudsol = SudokuSolver()
    s=sudsol.solve(Pole)
    if Pole == []:
        print('Nema riesenie')
    else:
        for n,r,c in s:
            print('{} {} {}'.format(n,r,c))
    
##    N = int(input())
##    nq = SudokuSolver()
##    s = nq.solve(N)
##    if len(s) == 0:
##        print('Nema riesenie')
##    else:
##        for r,c in s:
##            print('{} {}'.format(r,c))

# vim: set sw=4 ts=4 sts=4 et :
