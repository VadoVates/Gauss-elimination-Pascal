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

function EliminacjaGaussa (var macierz : tMacierz; ileNiewiadomych : word; eps : double; wektorX : tWektor) : boolean;
var
  i, j, kolumna, k, l : integer;
  mnoznik, x2 : double;
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
      for kolumna:=i+1 to ileNiewiadomych do
      macierz[j][kolumna] := macierz [j][kolumna] + mnoznik * macierz[i][kolumna];
    end;
  end;



    for i:=(ileNiewiadomych-1) downto 0 do
    begin
      x2:=macierz[i][ileNiewiadomych];
      for j:=ileNiewiadomych-1 downto (i+1) do
      begin
        x2:=x2-macierz[i][j] * wektorX[j];
        //jezeli przekatna ma ktorykolwiek element = 0, wtedy macierz jest osobliwa,
        //zwracamy false i konczymy funkcje
        if (Modul(macierz[i][i]) < eps) then
        begin
          EliminacjaGaussa := false;
          Exit;
        end;
        wektorX[i]:=x2/macierz[i][i];
      end;
      EliminacjaGaussa := true;
    end;

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

  write ('Podaj jakiej oczekujesz dokladnosci (bliskosc do zera): ');
  read (eps);

  //+1 kolumn, poniewaz dopisujemy do macierzy kolumne wyrazow wolnych
  SetLength (macierzAB, ileNiewiadomych, ileNiewiadomych+1);
  SetLength (wektorX,ileNiewiadomych);

{  //wypisywanie zeby sprawdzic czy dziala prawidlowo wczytanie i przetwarzanie macierzy
  for i:=0 to ileNiewiadomych-1 do
  begin
    //for j:=0 to ileNiewiadomych do
    write (X[i]:3:0);
    writeln;
  end;
}

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
    for i:=0 to ileNiewiadomych-1 do
    begin
      for j:=0 to ileNiewiadomych do
      write (macierzAB[i][j]:3:0);
      writeln;
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
