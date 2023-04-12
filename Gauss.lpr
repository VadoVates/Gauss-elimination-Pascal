program Gauss;

type
  tMacierz = array of array of double;
  tWektorDouble = array of double;
  tWektorIndeks = array of word;

function Modul (x:double) : double;
begin
  if (x<0) then
  Modul:=x*(-1)
  else
  Modul:=x;
end;

function EliminacjaGaussa (var macierz : tMacierz; ileNiewiadomych : word; eps : double;
                          var wektorX : tWektorDouble) : boolean;
var
  i, j, kolumna, k, l : word;
  mnoznik, suma : double;
begin
  //petla for - tworzymy macierz trojkatna gorna
  //do -2, bo nie potrzebujemy zerowac parametru przy ostatniej niewiadomej.
  for i:=0 to ileNiewiadomych-2 do
  begin
    //od i+1, wiec schodzi wiersz nizej niz i, a indeks i tworzy "schodki"
    for j:=i+1 to ileNiewiadomych-1 do
    begin
      //jezeli przekatna ma ktorykolwiek element = 0, wtedy macierz jest osobliwa,
      //zwracamy false i konczymy funkcje
      if (Modul(macierz[i][i]) < eps) then
      begin
        EliminacjaGaussa := false;
        Exit;
      end;
      mnoznik := (macierz[j][i]/macierz[i][i])*(-1);
      //dodawanie wiersza przemnozonego przez mnoznik
      for kolumna:=i to ileNiewiadomych do
      macierz[j][kolumna] := macierz [j][kolumna] + mnoznik * macierz[i][kolumna];

      //kontrolne wypisanie
      writeln ('Macierz AB:');
      for k:=0 to ileNiewiadomych-1 do
      begin
        for l:=0 to ileNiewiadomych do
        write (macierz[k][l]:8:3,' ');
        writeln;
      end;
      //
    end;
  end;

  //obliczanie x i wpisywanie wartosci x do wektora X
  //od konca, bo macierz trojkatna gorna
  for i:=(ileNiewiadomych-1) downto 0 do
  begin
    //suma jest rowna wyrazowi wolnemu
    suma := macierz[i][ileNiewiadomych];
    //j do i+1 po to, zeby petla nie wykonywala sie gdy wektorX jest pusty
    for j:=(ileNiewiadomych-1) downto (i+1) do
    begin
      suma := suma - macierz[i][j] * wektorX[j];
    end;
    //jezeli przekatna ma ktorykolwiek element = 0, wtedy macierz jest osobliwa,
    //zwracamy false i konczymy funkcje
    if (Modul(macierz[i][i]) < eps) then
    begin
      EliminacjaGaussa := false;
      Exit;
    end;
    wektorX[i]:=suma/macierz[i][i];
  end;
  EliminacjaGaussa := true;
end;

procedure Gauss (var macierzAB : tMacierz; ileNiewiadomych: word; eps : double);
var
  wektorX : tWektorDouble;
  i : word;
begin
  SetLength (wektorX,ileNiewiadomych);

  //jezeli funkcja zwrocila 'true', to oznacza, ze det!=0 i mozna wypisac wyniki
  if (EliminacjaGaussa (macierzAB, ileNiewiadomych, eps, wektorX)) then
  begin
    writeln ('Funkcja eliminacji Gaussa zwrocila true');
    for i:=0 to ileNiewiadomych-1 do
    begin
      writeln ('x', i+1, ' = ', wektorX[i]:8:4);
    end;
  end
  //jezeli funkcja zwrocila 'false', to oznacza, ze det=0 i nie mozna wypisac wynikow
  else
  begin
    writeln ('Macierz osobliwa, det = 0');
  end;
end;

function EliminacjaGaussaCrouta (var macierz : tMacierz; ileNiewiadomych : word; eps : double;
                          var wektorX : tWektorDouble; var wektorKolumna : tWektorIndeks) : boolean;
var
  i, j, kolumna, l : word;
  mnoznik, suma : double;
begin
  //petla for - tworzymy macierz trojkatna gorna
  //do -2, bo nie potrzebujemy zerowac parametru przy ostatniej niewiadomej.
  for i:=0 to ileNiewiadomych-2 do
  begin
    kolumna := i;
    //przeszukanie forem za liczba najdalsza od zera i zamiana miejscami kolumn
    for j:=i+1 to ileNiewiadomych-1 do
    begin
      if (Modul(macierz[i][wektorKolumna[kolumna]]) < Modul(macierz[i][wektorKolumna[j]])) then
         kolumna := j;
    end;
    j:=wektorKolumna[kolumna];
    wektorKolumna[kolumna]:=wektorKolumna[i];
    wektorKolumna[i]:=j;
    //sprawdzenie czy na przekatnej dalej jest jakies zero, zeby nie dzielic przez zero
    for j:=i+1 to ileNiewiadomych-1 do
    begin
      if (Modul(macierz[i][wektorKolumna[i]]) < eps) then
      begin
        EliminacjaGaussaCrouta := false;
        Exit;
      end;
      //wyliczenie mnoznika
      mnoznik:= (-1)*macierz[j][wektorKolumna[i]] / macierz [i][wektorKolumna[i]];

      //dodanie wiersza pomnozonego przez mnoznik do wiersza, na ktorym operujemy
      for kolumna:=i to ileNiewiadomych do
      begin
        macierz[j][wektorKolumna[kolumna]]:=macierz[j][wektorKolumna[kolumna]] + mnoznik * macierz[i][wektorKolumna[kolumna]];
      end;

      //kontrolne wypisanie
      writeln ('Macierz AB:');
      for l:=0 to ileNiewiadomych-1 do
      begin
        for kolumna:=0 to ileNiewiadomych do
        write (macierz[l][wektorKolumna[kolumna]]:8:3,' ');
        writeln;
      end;
      //
    end;
  end;

  //niewiadome
  for i:=ileNiewiadomych-1 downto 0 do
  begin
    if (Modul(macierz[i][wektorKolumna[i]])<eps) then
    begin
      EliminacjaGaussaCrouta:=false;
      Exit;
    end;
    suma:= macierz [i][ileNiewiadomych];
    for j:=ileNiewiadomych-1 downto i+1 do
    begin
      suma:= suma - macierz[i][wektorKolumna[j]]*wektorX[wektorKolumna[j]];
    end;
    wektorX[wektorKolumna[i]]:=suma/macierz[i][wektorKolumna[i]];
  end;
  EliminacjaGaussaCrouta := true;
end;

procedure GaussCrout (var macierzAB : tMacierz; ileNiewiadomych: word; eps : double);
var
  wektorX : tWektorDouble;
  wektorKolumna : tWektorIndeks;
  i : word;
begin
  SetLength (wektorX,ileNiewiadomych);
  SetLength (wektorKolumna, ileNiewiadomych+1);

  for i:=0 to ileNiewiadomych do
      wektorKolumna[i]:=i;

  //jezeli funkcja zwrocila 'true', to oznacza, ze det!=0 i mozna wypisac wyniki
  if (EliminacjaGaussaCrouta (macierzAB, ileNiewiadomych, eps, wektorX, wektorKolumna)) then
  begin
    writeln ('Funkcja eliminacji Gaussa-Crouta zwrocila true:');
    for i:=0 to ileNiewiadomych-1 do
    begin
      writeln ('x', i+1, ' = ', wektorX[i]:8:4);
    end;
  end
  //jezeli funkcja zwrocila 'false', to oznacza, ze det=0 i nie mozna wypisac wynikow
  else
  begin
    writeln ('Macierz osobliwa, det = 0');
  end;
end;

function EliminacjaGaussaWiersze (var macierz : tMacierz; ileNiewiadomych : word; eps : double;
                          var wektorX : tWektorDouble; var wektorWiersz : tWektorIndeks) : boolean;
var
  i, j, wiersz, l : word;
  mnoznik, suma : double;
begin
  for i:=0 to ileNiewiadomych-2 do
  begin
    wiersz := i;
    //przeszukanie forem za liczba najdalsza od zera i zamiana miejscami wierszy
    for j:=i+1 to ileNiewiadomych-1 do
    begin
      if (Modul(macierz[wektorWiersz[wiersz]][i]) < Modul(macierz[wektorWiersz[j]][i])) then
         wiersz := j;
    end;
    j:=wektorWiersz[wiersz];
    wektorWiersz[wiersz]:=wektorWiersz[i];
    wektorWiersz[i]:=j;
    //sprawdzenie czy na przekatnej dalej jest jakies zero, zeby nie dzielic przez zero
    for j:=i+1 to ileNiewiadomych-1 do
    begin
      if (Modul(macierz[wektorWiersz[i]][i]) < eps) then
      begin
        EliminacjaGaussaWiersze := false;
        Exit;
      end;
      //wyliczenie mnoznika
      mnoznik:= (-1)*macierz[wektorWiersz[j]][i] / macierz [wektorWiersz[i]][i];

      //dodanie wiersza pomnozonego przez mnoznik do wiersza, na ktorym operujemy
      for wiersz:=i to ileNiewiadomych do
      begin
        macierz[wektorWiersz[j]][wiersz]:=macierz[wektorWiersz[j]][wiersz] + mnoznik * macierz[wektorWiersz[i]][wiersz];
      end;

      //kontrolne wypisanie
      writeln ('Macierz AB:');
      for wiersz:=0 to ileNiewiadomych-1 do
      begin
        for l:=0 to ileNiewiadomych do
        write (macierz[wektorWiersz[wiersz]][l]:8:3,' ');
        writeln;
      end;
      //
    end;
  end;

  //obliczanie x i wpisywanie wartosci x do wektora X
  //od konca, bo macierz trojkatna gorna
  for i:=(ileNiewiadomych-1) downto 0 do
  begin
    //suma jest rowna wyrazowi wolnemu
    suma := macierz[wektorWiersz[i]][ileNiewiadomych];
    //j do i+1 po to, zeby petla nie wykonywala sie gdy wektorX jest pusty
    for j:=(ileNiewiadomych-1) downto (i+1) do
    begin
      suma := suma - macierz[wektorWiersz[i]][j] * wektorX[j];
    end;
    //jezeli przekatna ma ktorykolwiek element = 0, wtedy macierz jest osobliwa,
    //zwracamy false i konczymy funkcje
    if (Modul(macierz[wektorWiersz[i]][i]) < eps) then
    begin
      EliminacjaGaussaWiersze := false;
      Exit;
    end;
    wektorX[i]:=suma/macierz[wektorWiersz[i]][i];
  end;
  EliminacjaGaussaWiersze:=true;
end;

procedure GaussWiersze (var macierzAB : tMacierz; ileNiewiadomych : word; eps : double);
var
  wektorX : tWektorDouble;
  wektorWiersz : tWektorIndeks;
  i : word;
begin
  SetLength (wektorWiersz,ileNiewiadomych);
  SetLength (wektorX,ileNiewiadomych);

  for i:=0 to ileNiewiadomych-1 do
      wektorWiersz[i]:=i;

  //jezeli funkcja zwrocila 'true', to oznacza, ze det!=0 i mozna wypisac wyniki
  if (EliminacjaGaussaWiersze (macierzAB, ileNiewiadomych, eps, wektorX, wektorWiersz)) then
  begin
    writeln ('Funkcja eliminacji Gaussa z zamiana wierszy zwrocila true:');
    for i:=0 to ileNiewiadomych-1 do
    begin
      writeln ('x', i+1, ' = ', wektorX[i]:8:4);
    end;
  end
  //jezeli funkcja zwrocila 'false', to oznacza, ze det=0 i nie mozna wypisac wynikow
  else
  begin
    writeln ('Macierz osobliwa, det = 0');
  end;
end;

function EliminacjaGaussaJordana (var macierz : tMacierz; ileNiewiadomych : word; eps : double;
                          var wektorX : tWektorDouble; var wektorWiersz : tWektorIndeks) : boolean;
var
  i, j, wiersz, l : word;
  mnoznik, suma : double;
begin
  for i:=0 to ileNiewiadomych-2 do
  begin
    wiersz := i;
    //przeszukanie forem za liczba najdalsza od zera i zamiana miejscami wierszy
    for j:=i+1 to ileNiewiadomych-1 do
    begin
      if (Modul(macierz[wektorWiersz[wiersz]][i]) < Modul(macierz[wektorWiersz[j]][i])) then
         wiersz := j;
    end;
    j:=wektorWiersz[wiersz];
    wektorWiersz[wiersz]:=wektorWiersz[i];
    wektorWiersz[i]:=j;
    //sprawdzenie czy na przekatnej dalej jest jakies zero, zeby nie dzielic przez zero
    for j:=i+1 to ileNiewiadomych-1 do
    begin
      if (Modul(macierz[wektorWiersz[i]][i]) < eps) then
      begin
        EliminacjaGaussaJordana := false;
        Exit;
      end;
      //wyliczenie mnoznika
      mnoznik:= (-1)*macierz[wektorWiersz[j]][i] / macierz [wektorWiersz[i]][i];

      //dodanie wiersza pomnozonego przez mnoznik do wiersza, na ktorym operujemy
      for wiersz:=i to ileNiewiadomych do
      begin
        macierz[wektorWiersz[j]][wiersz]:=macierz[wektorWiersz[j]][wiersz] + mnoznik * macierz[wektorWiersz[i]][wiersz];
      end;

      //kontrolne wypisanie
      writeln ('Macierz AB:');
      for wiersz:=0 to ileNiewiadomych-1 do
      begin
        for l:=0 to ileNiewiadomych do
        write (macierz[wektorWiersz[wiersz]][l]:8:3,' ');
        writeln;
      end;
      //
    end;
  end;
  
  for i:=ileNiewiadomych-1 downto 1 do
  begin
    for j:=i-1 downto 0 do
    begin
      mnoznik:= (-1)*macierz[wektorWiersz[j]][i] / macierz [wektorWiersz[i]][i];
      for wiersz:=ileNiewiadomych downto j do
      begin
	macierz[wektorWiersz[j]][wiersz]:=macierz[wektorWiersz[j]][wiersz] + mnoznik * macierz[wektorWiersz[i]][wiersz];
      end;
    end;
    //kontrolne wypisanie
    writeln ('Macierz AB:');
    for wiersz:=0 to ileNiewiadomych-1 do
    begin
      for l:=0 to ileNiewiadomych do
      write (macierz[wektorWiersz[wiersz]][l]:8:3,' ');
      writeln;
    end;
    //
  end;

  //Wpisanie wynikow do wektora X wg kolejnosci wierszy i zrobienie jednostkowej macierzy
  for i:=ileNiewiadomych-1 downto 0 do
  begin
    macierz[wektorWiersz[i]][ileNiewiadomych]:= macierz[wektorWiersz[i]][ileNiewiadomych] / macierz [wektorWiersz[i]][i];
    macierz[wektorWiersz[i]][i] := macierz [wektorWiersz[i]][i] / macierz [wektorWiersz[i]][i];
    wektorX[i] := macierz [wektorWiersz[i]][ileNiewiadomych];
    //kontrolne wypisanie
    writeln ('Macierz AB:');
    for wiersz:=0 to ileNiewiadomych-1 do
    begin
      for l:=0 to ileNiewiadomych do
      write (macierz[wektorWiersz[wiersz]][l]:8:3,' ');
      writeln;
    end;
    //
  end;
  EliminacjaGaussaJordana:=true;
end;

procedure GaussJordan (var macierzAB : tMacierz; ileNiewiadomych : word; eps : double);
var
  wektorX : tWektorDouble;
  wektorWiersz : tWektorIndeks;
  i : word;
begin
  SetLength (wektorWiersz,ileNiewiadomych);
  SetLength (wektorX,ileNiewiadomych);
  for i:=0 to ileNiewiadomych-1 do
      wektorWiersz[i]:=i;

  //jezeli funkcja zwrocila 'true', to oznacza, ze det!=0 i mozna wypisac wyniki
  if (EliminacjaGaussaJordana (macierzAB, ileNiewiadomych, eps, wektorX, wektorWiersz)) then
  begin
    writeln ('Funkcja eliminacji Gaussa z zamiana wierszy zwrocila true:');
    for i:=0 to ileNiewiadomych-1 do
    begin
      writeln ('x', i+1, ' = ', wektorX[i]:8:4);
    end;
  end
  //jezeli funkcja zwrocila 'false', to oznacza, ze det=0 i nie mozna wypisac wynikow
  else
  begin
    writeln ('Macierz osobliwa, det = 0');
  end;
end;

procedure Jacobi (var macierz : tMacierz; ileNiewiadomych : word; eps : double);
var
  wektorX, wektorY : tWektorDouble;
  wektorWiersz : tWektorIndeks;
  i, wiersz, j, l : word;
  maxIter : word;
  suma, modulSuma : double;
  warunek : boolean;
begin
  SetLength (wektorWiersz,ileNiewiadomych);
  SetLength (wektorX,ileNiewiadomych);
  SetLength (wektorY,ileNiewiadomych);
  for i:=0 to ileNiewiadomych-1 do
      wektorWiersz[i]:=i;
  for i:=0 to ileNiewiadomych-1 do
  begin
    wiersz := i;
    //przeszukanie forem za liczba najdalsza od zera i zamiana miejscami wierszy
    for j:=i+1 to ileNiewiadomych-1 do
    begin
      if (Modul(macierz[wektorWiersz[wiersz]][i]) < Modul(macierz[wektorWiersz[j]][i])) then
         wiersz := j;
    end;
    j:=wektorWiersz[wiersz];
    wektorWiersz[wiersz]:=wektorWiersz[i];
    wektorWiersz[i]:=j;
    //sprawdzenie czy na przekatnej dalej jest jakies zero, zeby nie dzielic przez zero
    for j:=i+1 to ileNiewiadomych-1 do
    begin
      if (Modul(macierz[wektorWiersz[i]][i]) < eps) then
      begin
        writeln ('DET=0, nie da rady psze pana!');
      end;
    end;
    {kontrolne wypisanie
    writeln ('Macierz AB:');
    for wiersz:=0 to ileNiewiadomych-1 do
    begin
      for l:=0 to ileNiewiadomych do
      write (macierz[wektorWiersz[wiersz]][l]:8:3,' ');
      writeln;
    end;
    }
  end;

  for i:=0 to ileNiewiadomych-1 do
  begin
    macierz[wektorWiersz[i]][i] := 1/macierz[wektorWiersz[i]][i];
    macierz[wektorWiersz[i]][ileNiewiadomych] := macierz[wektorWiersz[i]][ileNiewiadomych] * macierz[wektorWiersz[i]][i];
    wektorY[i]:= macierz[wektorWiersz[i]][ileNiewiadomych];
  end;
  j:=0;
  for i:=0 to ileNiewiadomych-1 do
  begin
    for j:=0 to ileNiewiadomych-1 do
    begin
      if (j<>i) then macierz[wektorWiersz[i]][j]:=-macierz[wektorWiersz[i]][j]*macierz[wektorWiersz[i]][i];
    end;
  end;
  //kontrolne wypisanie
  writeln ('Macierz AB:');
  for wiersz:=0 to ileNiewiadomych-1 do
  begin
    for l:=0 to ileNiewiadomych do
    write (macierz[wektorWiersz[wiersz]][l]:8:3,' ');
    writeln;
  end;
  //
  writeln ('Dobra, kotles, ile iteracji?');
  read (maxIter);
  l:=0;

  repeat
    for i:=0 to ileNiewiadomych-1 do
    begin
      suma:=0;
      modulSuma:=0;
      for j:=0 to ileNiewiadomych-1 do
      begin
        if (j<>i) then
        begin
           suma:= suma + wektorY[j]*macierz[wektorWiersz[i]][j];
           modulSuma:= modulSuma + Modul(wektorY[j]*macierz[wektorWiersz[i]][j]);
        end;
        //if (Modul(macierz[wektorWiersz[i]][i])>modulSuma) then
        //   warunek:=true;
      end;
      wektorX[i]:=macierz[wektorWiersz[i]][ileNiewiadomych]+suma;
      writeln ('x',i+1,'=',wektorX[i]:8:4);
    end;
    for i:=0 to ileNiewiadomych-1 do
    begin
      wektorY[i]:=wektorX[i];
    end;
    l:=l+1;

  until (maxIter=l);
end;

procedure CzytajDane (var macierzAB : tMacierz; var ileNiewiadomych : word; var eps : double);
var
  i, j : word;
begin
  write ('Podaj ile chcesz niewiadomych: ');
  read (ileNiewiadomych);

  //write ('Podaj jakiej oczekujesz dokladnosci (bliskosc do zera): ');
  //read (eps);
  eps:=1e-12;
  //+1 kolumn, poniewaz dopisujemy do macierzy kolumne wyrazow wolnych
  //SetLength oprocz ustawienia wymiarow, wypelnia wszystkie komorki tabeli zerami
  SetLength (macierzAB, ileNiewiadomych, ileNiewiadomych+1);
  //numerowanie indeksow jest od 0, wiec liczymy do ile-1
  //odczytywanie zawartosci macierzy, tzn. parametrow
  for i:=0 to ileNiewiadomych-1 do
  begin
    for j:=0 to ileNiewiadomych do
    read (macierzAB[i][j]);
  end;
end;

procedure GotoweDane (var macierzAB : tMacierz; var ileNiewiadomych : word; var eps : double);
var
  i : word;
begin
  write ('Podaj czy wolisz dane latwiejsze [1], czy trudniejsze [2] (z zerami na przekatnej), czy [3]: ');
  read (i);
  ileNiewiadomych:=4;
  eps:=1e-12;
  SetLength (macierzAB, 4, 5);
  if (i=1) then
  begin
    macierzAB[0][0]:=4;
    macierzAB[0][1]:=-2;
    macierzAB[0][2]:=4;
    macierzAB[0][3]:=-2;
    macierzAB[0][4]:=8;
    macierzAB[1][0]:=3;
    macierzAB[1][1]:=1;
    macierzAB[1][2]:=4;
    macierzAB[1][3]:=2;
    macierzAB[1][4]:=7;
    macierzAB[2][0]:=2;
    macierzAB[2][1]:=4;
    macierzAB[2][2]:=2;
    macierzAB[2][3]:=1;
    macierzAB[2][4]:=10;
    macierzAB[3][0]:=2;
    macierzAB[3][1]:=-2;
    macierzAB[3][2]:=4;
    macierzAB[3][3]:=2;
    macierzAB[3][4]:=2;
  end
  else
  begin
    if (i=2) then
    begin
      macierzAB[0][0]:=0;
      macierzAB[0][1]:=2;
      macierzAB[0][2]:=3;
      macierzAB[0][3]:=4;
      macierzAB[0][4]:=49;
      macierzAB[1][0]:=1;
      macierzAB[1][1]:=0;
      macierzAB[1][2]:=3;
      macierzAB[1][3]:=4;
      macierzAB[1][4]:=45;
      macierzAB[2][0]:=1;
      macierzAB[2][1]:=2;
      macierzAB[2][2]:=0;
      macierzAB[2][3]:=4;
      macierzAB[2][4]:=36;
      macierzAB[3][0]:=1;
      macierzAB[3][1]:=2;
      macierzAB[3][2]:=3;
      macierzAB[3][3]:=0;
      macierzAB[3][4]:=23;
    end
    else
    begin
      if (i=3) then
      begin
        macierzAB[0][0]:=4;
        macierzAB[0][1]:=-1;
        macierzAB[0][2]:=-0.2;
        macierzAB[0][3]:=2;
        macierzAB[0][4]:=30;
        macierzAB[1][0]:=-1;
        macierzAB[1][1]:=5;
        macierzAB[1][2]:=0;
        macierzAB[1][3]:=-2;
        macierzAB[1][4]:=0;
        macierzAB[2][0]:=0.2;
        macierzAB[2][1]:=1;
        macierzAB[2][2]:=10;
        macierzAB[2][3]:=-1;
        macierzAB[2][4]:=-10;
        macierzAB[3][0]:=0;
        macierzAB[3][1]:=-2;
        macierzAB[3][2]:=-1;
        macierzAB[3][3]:=4;
        macierzAB[3][4]:=5;
      end
      else
      begin
        writeln ('Wyjscie');
        Exit;
      end;
    end;
  end;
end;

procedure Implementacja;
var
  macierzAB : tMacierz;
  ileNiewiadomych : word;
  eps : double;
  wybor : word;
begin
  repeat
    writeln ('Wybor metody');
    writeln ('1 - metoda Gaussa bazowa (problem z zerami na przekatnej)');
    writeln ('2 - metoda Gaussa z zamiana kolumn (Gaussa-Crouta)');
    writeln ('3 - metoda Gaussa z zamiana wierszy');
    writeln ('4 - metoda Gaussa-Jordana');
    writeln ('5 - metoda Jacobiego (iteracyjna)');
    writeln ('0 - wyjscie');
    readln (wybor);
    case (wybor) of
      1 :
      begin
        //CzytajDane (macierzAB, ileNiewiadomych, eps);
        GotoweDane (macierzAB, ileNiewiadomych, eps);
        writeln ('Metoda eliminacji Gaussa:');
        Gauss (macierzAB, ileNiewiadomych, eps);
      end;
      2 :
      begin
        //CzytajDane (macierzAB, ileNiewiadomych, eps);
        GotoweDane (macierzAB, ileNiewiadomych, eps);
        writeln ('Metoda eliminacji Gaussa-Crouta:');
        GaussCrout (macierzAB, ileNiewiadomych, eps);
      end;
      3 :
      begin
        //CzytajDane (macierzAB, ileNiewiadomych, eps);
        GotoweDane (macierzAB, ileNiewiadomych, eps);
        writeln ('Metoda eliminacji Gaussa z zamiana wierszy:');
        GaussWiersze (macierzAB, ileNiewiadomych, eps);
      end;
      4 :
      begin
        //CzytajDane (macierzAB, ileNiewiadomych, eps);
        GotoweDane (macierzAB, ileNiewiadomych, eps);
        writeln ('Metoda eliminacji Gaussa-Jordana:');
        GaussJordan (macierzAB, ileNiewiadomych, eps);
      end;
      5 :
      begin
        //CzytajDane (macierzAB, ileNiewiadomych, eps);
        GotoweDane (macierzAB, ileNiewiadomych, eps);
        writeln ('Metoda eliminacji Jacobiego (iteracyjna):');
        Jacobi (macierzAB, ileNiewiadomych, eps);
      end;
      0 :
      begin
        writeln ('Wyjscie z programu');
        Exit;
      end;
      else
      begin
        writeln ('Wpisales niepoprawny znak');
      end;
    end;
  until (wybor=0);
end;

begin
  Implementacja;
  writeln ('Wcisnij cokolwiek zeby kontynuowac wyjscie z programu');
  readln;
end.
