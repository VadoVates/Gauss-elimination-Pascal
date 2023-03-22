program Gauss;

type
  tMacierz = array of array of double;
  tWektor = array of double;

function Modul (x:double) : double;
begin
  if (x<0) then
  Modul:=x*(-1)
  else
  Modul:=x;
end;

function EliminacjaGaussa (var macierz : tMacierz; ileNiewiadomych : word;
                          eps : double; var wektorX : tWektor) : boolean;
var
  a , b : integer;
  i, j, kolumna : word;
  mnoznik, suma : double;
begin
  for i:=0 to ileNiewiadomych-2 do
  begin
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
	  
	  //
	  //
	  writeln ('Macierz AB:');
      for a := 0 to ileNiewiadomych - 1 do
      begin
        for b := 0 to ileNiewiadomych  do write ( macierz [ a ][ b ]:8:3 );
        writeln;
      end;
	  //
	  //
	  
    end;
  end;

  for i:=(ileNiewiadomych-1) downto 0 do
  begin
    suma := macierz[i][ileNiewiadomych];
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
	
	//
	//
	writeln ('Wektor X:');
    for a := 0 to ileNiewiadomych - 1 do
    begin
      write ( wektorX [ a ]:8:3 );
      writeln;
    end;
	//
	//
	
  end;
  EliminacjaGaussa := true;
end;

procedure Implementacja;
var
  macierzAB : tMacierz;
  wektorX : tWektor;
  ileNiewiadomych, i, j : word;
  eps : double;
begin
  write ('Podaj ile chcesz niewiadomych: ');
  read (ileNiewiadomych);

  //write ('Podaj jakiej oczekujesz dokladnosci (bliskosc do zera): ');
  //read (eps);
  eps:=1e-12;
  //+1 kolumn, poniewaz dopisujemy do macierzy kolumne wyrazow wolnych
  //SetLength oprocz ustawienia wymiarow, wypelnia wszystkie komorki tabeli zerami
  SetLength (macierzAB, ileNiewiadomych, ileNiewiadomych+1);
  SetLength (wektorX,ileNiewiadomych);

  //numerowanie indeksow jest od 0, wiec liczymy do ile-1
  //odczytywanie zawartosci macierzy, tzn. parametrow
  for i:=0 to ileNiewiadomych-1 do
  begin
    for j:=0 to ileNiewiadomych do
    read (macierzAB[i][j]);
  end;

  //jezeli funkcja zwrocila 'true', to oznacza, ze det!=0 i mozna wypisac wyniki
  if (EliminacjaGaussa (macierzAB, ileNiewiadomych, eps, wektorX)) then
  begin
    writeln ('Funkcja zwrocila trv');
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

begin
  Implementacja;
  readln;
  readln;
end.
