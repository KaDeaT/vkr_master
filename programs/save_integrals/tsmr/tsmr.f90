module tsmr
implicit none
contains
subroutine tsmr_integ(n,u,na,np,ln,rtol,atol,fpar, k_planets)
!Глобальные константы:  vtn - число затабулированных значений для v^{-1}, u^{-1}
!                       pmin,pmax,pdif - минимальный, максимальный и разница порядков для оценки остатка
integer,parameter :: vtn=99, pmin=5,pdif=2,pmax=60+pdif,fid=307
real,parameter :: dhi=5.0,dhd=0.2 !границы изменения длины шага для переключения порядка
!Глобальные переменные: n - размерность системы, u - число мономов, na - число ненулевых коэффициентов
!                       np - число точек вывода решения, sch - схема, ia,ja - индексы ненулевых коэффициентов,
!                       deg - степень правой части
!Локальные переменные: pwork - текущий порядок, i,j,ix - вспомогательные, ipar_alp - целочисленные параметры для alpha
integer :: n,u,na,np,sch(2*u),ia(n+1),ja(na),deg,pwork,i,j,ipar_alp(1),ix, k_planets
!Глобальные переменные: a,aa - ненулевые коэффициенты, их модули, vtb - значения v^{-1} соответствующие степени правой части,
!                       rdeg=1/deg, atol,rtol - требуемые абсолютная и относительная погрешности, af - вектор свободный членов,
!                       ut - вспомогательный массив
!Локальные переменные:  tx - решение, tx1,tx2 - вспомогательные массивы, x - коэффициенты тейлора, ox - точки вывода решения,
!                       time - времена расчета коэффициетов тейлора, tc,te - текущая и конечная точки интегрирования,
!                       h,ht - величина шага, величина шага при смене порядка, alp - масш. множители, dir - направление интегрирования,
!                       rpar_alp - вещественные параметры для alpha 
real :: tx(n),tx1(n),tx2(n),x((u+n+1)*(pmax+1)),a(na),aa(na),ox(np),vtb(vtn*(pmax-pmin)),rdeg,atol,rtol,&
        af(n),time(pmax-(pmin+pdif)+1),ut(n+u+1),tc,te,h,ht,alp(n),dir,rpar_alp(1),ho,& 
		mass(10), h_0, h_1, h_error, x0(n), h_c0(3), h_c1(3)
!Локальные переменные: fpar - имена файлов с данными        
character(*) :: fpar(:)
!Локальные переменные: ln - линейность системы, last - для выдачи, последний шаг
logical :: ln,last
!чтение данных
call tsmr_data(tx,mass,fpar)

x0 = tx

x=0.0
x(1)=1.0
aa=abs(a)
tc=ox(1)
te=ox(np)
dir=sign(1.0,te-tc)
last=.false.
!замер времени
call tsmr_time(x,tx)
!вычисление начального приближения к шагу
pwork=pmin+pdif
call tsmr_ualp(alp,tx,rpar_alp,ipar_alp)
h=tsmr_guar(tx,alp,pwork)*dir
!вычисление оптимального порядка и первого шага
call tsmr_order(pwork,h,tx,x,alp,.true.,tx1,tx2)

ht=h
!файл для записи
open(unit=fid,file=trim(fpar(5)))
write(fid,*) ox(1) 
write(fid,'(f36.30)') (tx(i),i=1,n)

open(unit=1992, file="proverka.txt", status="replace")
56 format(f36.30)
57 format(e37.30)

h_0 = calc_he(tx, n, mass, k_planets)
call calc_hc(tx, h_c0,mass, k_planets)
write(1992,*) 'h_energy0'
write(1992,56) h_0 
write(1992, *) 'h_square0'
write(1992,56) (h_c0(i), i=1,3)

ix=2
call tsmr_coef(x,tx,0,pmax)
do  
    !вычисление шага
    h=tsmr_step(h,tx,x,alp,pwork,tx1,tx2)       
    !проверка, вычисление порядка
    if(h/ht<dhd .or. h/ht>dhi) then
        call tsmr_order(pwork,h,tx,x,alp,.false.,tx1,tx2)
        ht=h
    end if   
    h=tsmr_step(h,tx,x,alp,pwork,tx1,tx2)        
    !выдача
    if((tc+h-ox(ix))*dir>0) then
        do 
            ho=ox(ix)-tc
            call tsmr_val(tx1,x,ho,0,pwork)
            write(fid, *) ox(ix)
			write(fid,'(f36.30)') (tx1(i),i=1,n) 
			
			write(1992,*)"coordinates"
			call rell_error(x0,tx1,n)

			h_1 = calc_he(tx1, n, mass, k_planets)
			write(1992,*) 'h_energy', ix
			write(1992,56) h_1 
			h_error = abs(h_1 - h_0)/abs(h_0)
			write(1992,*) 'h_energy error:'
			write(1992,57) h_error 
			
			call calc_hc(tx1,h_c1,mass,k_planets)
			write(1992, *) 'h_square'
			write(1992,56) (h_c1(i), i=1,3)
			write(1992,*) 'h_square error:'
			call rell_error(h_c0,h_c1,3)
			
            if(ix==np) then
                last=.true.
                exit
            else
                ix=ix+1
            end if
            if((tc+h-ox(ix))*dir<0) exit
        end do
    end if
    if(last) then
        exit
    else
        call tsmr_val(tx,x,h,0,pwork-pdif)
        call tsmr_coef(x,tx,0,pwork)
        tc=tc+h
    end if
end do
close(fid)
contains

subroutine tsmr_coef(x,tx,pl,pu)
! КОЭФФИЦИЕНТЫ ТЕЙЛОРА
! -------------ВХОДНЫЕ----------------------------
! x - массив для коэффициентов Тейлора
! tx - массив значений функции
! pl - порядок, от которого считаются коэффициенты/нужно с 0.
! pu - порядок, до которого считаются
!--------------ГЛОБАЛЬНЫЕ-------------------------
! n - размерность
! u - число мономов
! pmax - максимальный порядок
! a,ia,ja - CSR матрица коэффициентов системы
! sch - схема
!--------------ВОЗВРАЩАЕТ-------------------------
!x - коэффициенты Тейлора
integer :: pl,pu,i,j,k,m1,m2,m3,m4
real :: tx(n),x((u+n+1)*(pmax+1)),r,t
m1=pmax+1
do i=1,n
	x(m1*i+1)=tx(i)
end do
do i=pl,pu-1
	r=1.0/real(i+1)
	do j=1,u
		m2=2*j
		m3=sch(m2-1)*m1+1
		m4=sch(m2)*m1+1
		x((j+n)*m1+i+1)=dot_product(x(m3:m3+i),x(m4+i:m4:-1))
	end do
	do j=1,n
		t=0.0
		do k=ia(j),ia(j+1)-1
			t=t+a(k)*x(ja(k)*m1+i+1)
		end do
		x(j*m1+i+2)=r*t
	end do
end do	
end subroutine tsmr_coef

subroutine tsmr_val(tx,x,h,pl,pu)
! ЗНАЧЕНИЕ ПОЛИНОМОВ ТЕЙЛОРА
! -------------ВХОДНЫЕ----------------------------
! tx - массив для значения функции
! x - массив коэффициентов Тейлора
! h - длина шага
! pl - порядок, от которого считаются коэффициенты/нужно с 0.
! pu - порядок, до которого считаются
!--------------ГЛОБАЛЬНЫЕ-------------------------
! n - размерность
! u - число мономов
! pmax - максимальный порядок
!--------------ВОЗВРАЩАЕТ-------------------------
! tx - значение функции
integer::pl,pu,i,j,m
real::tx(n),x((u+n+1)*(pmax+1)),h
m=pmax+1
tx=0.0
do i=pu,pl,-1
	do j=1,n
		tx(j)=tx(j)*h+x(j*m+i+1)
	end do
end do
if(pl/=0) tx=tx*(h**pl)
end subroutine tsmr_val


real function tsmr_guar(tx,alp,pwork)
! ВЫЧИСЛЕНИЕ ШАГА ПО МАСШ. МНОЖИТЕЛЯМ
! -------------ВХОДНЫЕ----------------------------
! tx - массив значений функции
! alp - массив масш. множителей
! pwork - текущий порядок
!--------------ГЛОБАЛЬНЫЕ-------------------------
! n - размерность
real :: tx(n),alp(n),r,v
integer :: i,pwork
r=tsmr_rest(alp)
v=tsmr_inv(alp,tx,r,pwork-pdif)
tsmr_guar=v*r
end function tsmr_guar

real function tsmr_heur(h0,tx,x,pwork,tx1,tx2)
! ВЫЧИСЛЕНИЕ ШАГА СТАНДАРТНОЙ КОРРЕКЦИЕЙ
! -------------ВХОДНЫЕ----------------------------
! h0 - приближение к длине шага
! tx - массив значений функции
! x - массив коэффициентов Тейлора
! pwork - текущий порядок   
! tx1,tx2 - вспомогательные массивы
!--------------ГЛОБАЛЬНЫЕ-------------------------
! n - размерность
! u - число мономов
! pmax - максимальный порядок
! pdif - разница порядков
! atol, rtol - требуемые погрешности   
integer :: pwork
real :: tx(n),tx1(n),tx2(n),x((u+n+1)*(pmax+1)),h0,s,t
t=0.0
call tsmr_val(tx1,x,h0,0,pwork-pdif)
call tsmr_val(tx2,x,h0,pwork-pdif+1,pwork)
do i=1,n
	s=atol+max(abs(tx(i)),abs(tx1(i)))*rtol
	t=t+(tx2(i)/s)**2
end do
if(t==0.0) then
	tsmr_heur=h0
else
	t=sqrt(t/real(n))
	tsmr_heur=h0*(1.0/t)**(1.0/real(pwork-pdif+1))
end if
end function tsmr_heur


real function tsmr_corr(h0,x,pwork,tx1,tx2)
! ВЫЧИСЛЕНИЕ ШАГА ИТЕРАТИВНОЙ КОРРЕКЦИЕЙ
! -------------ВХОДНЫЕ----------------------------
! h0 - приближение к длине шага
! x - массив коэффициентов Тейлора
! pwork - текущий порядок   
! tx1,tx2 - вспомогательные массивы
!--------------ГЛОБАЛЬНЫЕ-------------------------
! n - размерность
! u - число мономов
! pmax - максимальный порядок
! atol, rtol - требуемые погрешности 
! pdif - разница порядков
integer,parameter :: d=3
real :: tx(n),tx1(n),tx2(n),x((u+n+1)*(pmax+1)),h0,h,dh,e,s1,s2,s3
integer :: pwork,i
dh=h0/real(d)
call tsmr_val(tx1,x,h0,0,pwork-pdif)
call tsmr_val(tx2,x,h0,pwork-pdif+1,pwork)
e=maxval(abs(tx2)/(abs(tx1)+atol)) 
s1=sign(1.0,rtol-e)
i=1
do 
	h=h0+s1*real(i)*dh
	call tsmr_val(tx1,x,h,0,pwork-pdif)
	call tsmr_val(tx2,x,h,pwork-pdif+1,pwork)
	e=maxval(abs(tx2)/(abs(tx1)+atol))

	s2=sign(1.0,rtol-e)

	if(s1*s2<0.0) then
        s3=real(i)
        if(rtol<=e) s3=real(i-1)
		tsmr_corr=h0+s1*s3*dh
		exit
	else if(i==d-1) then
		h0=h; dh=h0/real(d); i=1
		cycle
	end if
	i=i+1
end do
end function tsmr_corr

real function tsmr_step(h0,tx,x,alp,pwork,tx1,tx2)
! АВТОМАТИЧЕСКОЕ ВЫЧИСЛЕНИЕ ШАГА 
! -------------ВХОДНЫЕ----------------------------
! h0 - приближение к длине шага
! tx - массив значений функции
! x - массив коэффициентов Тейлора
! alp - массив для масш. множителей
! pwork - текущий порядок
! tx1,tx2 - вспомогательные массивы
!--------------ГЛОБАЛЬНЫЕ-------------------------
! n - размерность
! u - число мономов
! dir - направление интегрирования
! pmax - максимальный порядок
! rpar_alp, ipar_alp - параметры для масш. множителей
real::h0,tx(n),x((u+n+1)*(pmax+1)),tx1(n),tx2(n),alp(n),h1,h2,h3
integer :: pwork
!прибилижение
call tsmr_ualp(alp,tx,rpar_alp,ipar_alp)
h1=tsmr_guar(tx,alp,pwork) 
h2=tsmr_heur(h0,tx,x,pwork,tx1,tx2)
h3=max(abs(h1),abs(h2))*dir
!корректировка
tsmr_step=tsmr_corr(h3,x,pwork,tx1,tx2)
end function tsmr_step

subroutine tsmr_data(tx, mass, fpar)
! ЧТЕНИЕ ДАННЫХ ИЗ ФАЙЛОВ 
! -------------ВХОДНЫЕ----------------------------
! tx - массив для значений функции
! mass -- массив для масс планет
! fpar(1) - файл с начальными данными
! fpar(2) - файл со схемой
! fpar(3) - файл с коэффициентами
! fpar(4) - файл с таблицей
! fpar(6) - файл с точками выдачи
! fpar(7) - файл с массами планет
!--------------ГЛОБАЛЬНЫЕ-------------------------
! n - размерность
! u - число мономов
! na - число коэффициентов 
! np - число точек выдачи
! ln - линейность
! a,ia,ja - массивы с коэффициентами и индексами
! sch - массив со схемой
! ox - массив с точками выдачи
! rdeg - 1/(степень правой части)
! vtn - число точек на порядок в таблице
! vtb - таблица
!--------------ВОЗВРАЩАЕТ-------------------------
! tx - массив значений решения
! a,ia,ja - массивы с коэффициентами и индексами
! sch - массив со схемой
! ox - массив с точками выдачи
! vtb - таблица 
integer,parameter :: fid=301
integer :: i,k1,k2,l,it(na), k_planets 
real :: tx(:), t, mass(10)
character(*)::fpar(:)
! начальные данные
open(unit=fid,file=trim(adjustl(fpar(1))),status='old')
read(fid,*) (tx(i),i=1,n)
close(fid)

! массы планет
open(unit=fid, file=trim(adjustl(fpar(7))), status="old")
k_planets = int((-11+int(sqrt(real(169+8*n))))/2)
mass(k_planets) = 1.0d0
read(fid,*) (mass(i), i = 1, k_planets-1)
close(fid)

! ненулевые коэффициенты, свободные члены
open(unit=fid,file=trim(fpar(3)),status='old')
k1=-1
k2=0
read(fid,*) (a(i),it(i),ja(i),i=1,na)
ia(n+1)=na+1
do i=1,na
    l=it(i)
	if(ja(i)==0) then
		af(l)=a(i)
	end if
	if(l>k1) then
		k2=k2+1
		ia(k2)=i		
		k1=l
    end if 
end do
close(fid)

if(ln) then
  deg=1.0
	rdeg=1.0
else
    open(unit=fid,file=trim(fpar(2)),status='old')
    read(fid,*) (sch(2*i-1),sch(2*i),i=1,u) 
    close(fid)
    deg=tsmr_deg()
    rdeg=1.0/real(deg-1)
end if
! таблица

open(unit=fid,file=trim(fpar(4)),status='old')
do i=1,(deg-1)*vtn*100+(pmin-1)*vtn
	read(fid,*) 
end do
do i=1,vtn*(pmax-pmin) 
	read(fid,*) vtb(i)
end do
close(fid)

! точки вывода
open(unit=fid,file=trim(fpar(6)),status='old')
read(fid,*) (ox(i),i=1,np)
close(fid)
end subroutine tsmr_data

integer function tsmr_deg()
! СТЕПЕНЬ ПРАВОЙ ЧАСТИ 
!--------------ГЛОБАЛЬНЫЕ-------------------------
! sch - массив со схемой
integer :: i,t(n+u),s1,s2
t(1:n)=1
do i=1,u
	s1=sch(2*i-1)
	s2=sch(2*i)
	t(n+i)=t(s1)+t(s2)
end do
tsmr_deg=maxval(t)
end function tsmr_deg

real function tsmr_inv(alp,tx,r,pwork)
! ВЫЧИСЛЕНИЕ TAU
! -------------ВХОДНЫЕ----------------------------
! alp - массив масш. множителей
! tx - массив для значений функции
! r - оценка радиуса сходимости
! pwork  -текущий порядок
!--------------ГЛОБАЛЬНЫЕ-------------------------
! n - размерность
! rtol,atol - погрешности
! pmin,pmax - мин. и макс. порядки 
! vtn - число точек на порядок в таблице
! vtb - таблица
real,parameter :: v=0.01
real :: alp(n),tx(n),e,s,r
integer :: pwork,i,q
if(ln) then
	s=minval((abs(tx)+atol)/alp)
	e=rtol*s
	s=maxval((abs(tx)+atol)/alp)
	e=e/(s+r*maxval(af/alp))
else
	s=minval((abs(tx)+atol)/alp)
	e=rtol*s
end if
i=(pwork-pmin)*vtn
q=tsmr_vsearch(e,vtb(i+1:i+vtn))
if(q==1) then
	tsmr_inv=v*e/vtb(1)
else
	i=q-1
	tsmr_inv=v*(e+i*vtb(q)-q*vtb(i))/(vtb(q)-vtb(i))
end if
end function tsmr_inv

integer function tsmr_vsearch(e,v)
! ПОИСК В МАССИВЕ
! -------------ВХОДНЫЕ----------------------------
! e - число
! v - массив
integer::i,j,k
real::e,v(vtn)
i=1
j=vtn 
DO 
	k=(i+j)/2
	IF(e<v(k)) THEN 
		j=k  
	ELSE
		i=k
	END IF
	IF (i+1>=j) then
		k=i
		EXIT	
	end if
END DO
tsmr_vsearch=k
end function tsmr_vsearch

subroutine tsmr_time(x,tx)
! ЗАМЕР ВРЕМЕНИ
! -------------ВХОДНЫЕ----------------------------
! x - массив для коэффициентов Тейлора
! tx - массив значений функции
!--------------ГЛОБАЛЬНЫЕ-------------------------
! n - размерность
! u - число мономов
! pmin,pmax - мин. и макс. порядки 
! time - массив для значений времени
!--------------ВОЗВРАЩАЕТ-------------------------
! time - массив времен рассчета коэффициентов
real::x((u+n+1)*(pmax+1)),tx(n),tx1(n)
real(kind=4) :: t,ts,tf
real(kind=4),parameter :: dt=0.1
real,parameter :: h=0.001
integer :: i,j,k
real :: s
do i=1,pmax-pmin-pdif+1
	j=1
	k=i-1+pmin+pdif
	call cpu_time(ts)
	do	
		call tsmr_coef(x,tx,0,k)
		!~ call tsmr_val(tx1,x,h,0,k)
		call cpu_time(tf)
		t=tf-ts
		if(t>dt) then
      s=real(t)  
			time(i)=s/real(j+1)
			j=1
			exit
		else
			j=j+1
		end if
	end do	
end do
end subroutine tsmr_time

subroutine tsmr_order(pwork,h0,tx,x,alp,first,tx1,tx2)
! ВЫБОР ПОРЯДКА
! -------------ВХОДНЫЕ----------------------------
! pwork - текущий порядок
! h0 - текущая длина шага
! tx - массив значений функции
! x - массив для коэффициентов Тейлора
! alp - массив масш. множителей
! first - первый шаг
!--------------ГЛОБАЛЬНЫЕ-------------------------
! n - размерность
! u - число мономов
! pmax - макс. порядок
! time - массив времен рассчета коэффициентов
!--------------ВОЗВРАЩАЕТ-------------------------
! pwork - новый порядок
! h - соответствующая длина шага
real::v1,v2,h1,h2,h0,alp(n),tx(n),tx1(n),tx2(n),x((u+n+1)*(pmax+1))
integer::i,p,pwork
logical::first,cont
p=pwork
cont=.true.
if(first) then
    v2=0.0
    call tsmr_coef(x,tx,0,pmax)
    do i=pmin+pdif,pmax
        h1=tsmr_step(h0,tx,x,alp,i,tx1,tx2)
        v1=h1/time(i-pmin-pdif+1)
        if(v1>v2) then
            v2=v1
            h2=h1
            p=i
        end if
    end do
else
    v2=h/time(pwork-pmin-pdif+1)
    do i=pwork-1,pmin+pdif,-1
        h1=tsmr_step(h0,tx,x,alp,i,tx1,tx2)
        v1=h1/time(i-pmin-pdif+1)
        if(v1>v2) then
            h2=h1
            p=i
            cont=.false.
            exit
        end if
    end do
    if(cont) then
        do i=pwork+1,pmax
            call tsmr_coef(x,tx,i-1,i)
            h1=tsmr_step(h0,tx,x,alp,i,tx1,tx2)
            v1=h1/time(i-pmin-pdif+1)
            if(v1>v2) then
                h2=h1
                p=i
                exit
            end if
        end do
    end if
end if
pwork=p
h=h2
end subroutine tsmr_order

real function tsmr_rest(alp)
! ОЦЕНКА РАДИУСА СХОДИМОСТИ
! -------------ВХОДНЫЕ----------------------------
! alp - массив масш. множителей
!--------------ГЛОБАЛЬНЫЕ-------------------------
! n - размерность
! u - число мономов
! ut - вспомогательный массив
! sch - массив со схемой
! rdeg - 1/(степень правой части)
! aa,ia,ja - модули коэффициентов и их индексы
integer:: i,j,m
real::alp(n),s,s1
! вычисление правой части
ut(1)=1.0
ut(2:n+1)=alp
do i=1,u
	m=2*i
	ut(i+n+1)=ut(sch(m-1))*ut(sch(m))
end do
s1=0.0
do i=1,n
	s=0.0
	do j=ia(i),ia(i+1)-1
		s=s+aa(j)*ut(ja(j)+1)
	end do
	s=s/alp(i)
	if(s>s1) s1=s
end do
tsmr_rest=rdeg/s1
end function tsmr_rest


subroutine tsmr_ualp(alp,tx,rpar,ipar)
! МАСШТАБИРУЮЩИЕ МНОЖИТЕЛИ
! -------------ВХОДНЫЕ----------------------------
! alp - массив масш. множителей
! tx - массив значений функции
! rpar, ipar - доп. массивы с параметрами для alpha
!--------------ГЛОБАЛЬНЫЕ-------------------------
! n - размерность
real :: alp(n),tx(n),rpar(:),t
integer :: ipar(:)
optional :: rpar,ipar
t=maxval(abs(tx))
if(t==0.0) t=1.0
alp=1.0/t
end subroutine tsmr_ualp


subroutine rell_error(x,y,k)
! Относительные погрешности 
! -------------ВХОДНЫЕ----------------------------
! x - массив нулевых значений
! y - массив значений в текущей точке
!--------------ГЛОБАЛЬНЫЕ-------------------------
! k - размерность
integer :: k
real :: x(k), y(k), temp(k), max_error
57 format(e37.30)
temp = abs(x-y)/abs(x)
max_error = maxval(temp)
write(1992,*) "rell_error:"
write(1992,57) max_error
end subroutine rell_error

!------------------------------------------------------------------------------
! calc_he()	-- функция вычисляющая значение интеграла энергии в точке x
!------------------------------------------------------------------------------
! n 		-- размерность x
! k_planets	-- количество планет в задаче N-тел. Солнце тоже учитывается.
! x 		-- вектор х (координаты и скорости планет)
! mass		-- массы планет (относительно Солнца. Масса Солнца = 1) 
!------------------------------------------------------------------------------
! gamma		-- квадрат постоянной Гаусса (gamma = k^2, k = 0,01720209895) 
!------------------------------------------------------------------------------
function calc_he(x, n, mass, k_planets)
integer :: i, j, k, n, k_planets
real, parameter :: 	gamma = 0.0002959122082855911025
real :: x(n), mass(10), mass_obsh, r_p(45), r_s(3), T, U, S3, calc_he 

	r_s = 0
	r_p = 0	
	mass_obsh = 0
	
	do i = 1, k_planets
		mass_obsh = mass_obsh + mass(i)
	end do

! Вычисление r_{0i}	
	do i = 1, k_planets-1   
		r_p(i) = sqrt(x(3*(i-1)+1)**2 + x(3*(i-1)+2)**2 + x(3*(i-1)+3)**2)
	end do

! Вычисление r_{ij}
	k = 1
	do i = 1, k_planets-2
		do j = i+1, k_planets-1
			r_p(k_planets-1 + k) = sqrt((x(3*(i-1)+1) - x(3*(j-1)+1))**2 + (x(3*(i-1)+2) - x(3*(j-1)+2))**2 + (x(3*(i-1)+3) - x(3*(j-1)+3))**2) 
		end do
	end do

! Вычисление потенципльной составляющей
	do i = 1,k_planets-1   
		U = U + mass(i)/r_p(i)
	end do

	k = 1
	do i = 1, k_planets-2
		do j = i+1, k_planets-1
			U = U + (mass(i)*mass(j))/r_p(k_planets-1+k) 
		end do
	end do
	U = gamma*U

! Вычисление кинетической составляющей
	do i = 1,k_planets-1   
		T = T + mass(i)*(x(3*(k_planets + i) - 6 + 1)**2 + x(3*(k_planets + i) - 6 + 2)**2 + x(3*(k_planets + i) - 6 + 3)**2)
	end do
	T = 0.5*T

	do i = 1,k_planets-1   
		r_s(1) = r_s(1) + mass(i)*x(3*(k_planets + i) - 6 + 1)
		r_s(2) = r_s(2) + mass(i)*x(3*(k_planets + i) - 6 + 2)
		r_s(3) = r_s(3) + mass(i)*x(3*(k_planets + i) - 6 + 3)
	end do
	S3 = (1/(2*mass_obsh))*(r_s(1)**2 + r_s(2)**2 + r_s(3)**2)

	calc_he = T - U - S3
end function calc_he

subroutine calc_hc(x,h_c,mass, k_planets)
integer :: i, k_planets
real :: x(n), h_c(3), mass(10), mass_obsh, mx, my, mz, mxdot, mydot, mzdot 

	mass_obsh = 0
	h_c = 0	
	mx = 0
	my = 0
	mz = 0
	mxdot = 0
	mydot = 0
	mzdot = 0
	do i = 1, k_planets
		mass_obsh = mass_obsh + mass(i)
	end do
!	write(1992,*) "m_obsh", mass_obsh
	
	do i = 1,k_planets-1   
		h_c(1) = h_c(1) +  mass(i)*(x(3*(i-1)+2)*x(3*(k_planets + i) - 6 + 3) - x(3*(i-1)+3)*x(3*(k_planets + i) - 6 + 2)) 
		h_c(2) = h_c(2) +  mass(i)*(x(3*(i-1)+3)*x(3*(k_planets + i) - 6 + 1) - x(3*(i-1)+1)*x(3*(k_planets + i) - 6 + 3))
		h_c(3) = h_c(3) +  mass(i)*(x(3*(i-1)+1)*x(3*(k_planets + i) - 6 + 2) - x(3*(i-1)+2)*x(3*(k_planets + i) - 6 + 1))
		mx = mx + mass(i)*x(3*(i-1)+1)
		my = my + mass(i)*x(3*(i-1)+2)
		mz = mz + mass(i)*x(3*(i-1)+3)
		mxdot = mxdot + mass(i)*x(3*(k_planets + i) - 6 + 1)
		mydot = mydot + mass(i)*x(3*(k_planets + i) - 6 + 2)
		mzdot = mzdot + mass(i)*x(3*(k_planets + i) - 6 + 3)
	end do
	
!	write(1992,*) "hc", h_c(1), h_c(2), h_c(3)
!	write(1992,*) "mx, mxdot", mx, mxdot
!	write(1992,*) "my, mydot", my, mydot
!	write(1992,*) "mz, mzdot", mz, mzdot
	
	h_c(1) = h_c(1) - (1/mass_obsh)*(my*mzdot - mz*mydot)
	h_c(2) = h_c(2) - (1/mass_obsh)*(mz*mxdot - mx*mzdot)
	h_c(3) = h_c(3) - (1/mass_obsh)*(mx*mydot - my*mxdot)
end subroutine calc_hc

end subroutine tsmr_integ

subroutine tsmr_cfg(cfg,ipar,rpar,ln,fpar)
! ЧТЕНИЕ КОНФИГУРАЦИОННОГО ФАЙЛА
! -------------ВХОДНЫЕ----------------------------
! cfg - имя конфигурационного файла
! ipar - массив с целыми числами
! rpar - массив с вещественными числами
! ln - логическая переменная
! fpar(1) - файл с начальными данными
! fpar(2) - файл со схемой
! fpar(3) - файл с коэффициентами
! fpar(6) - файл с точками выдачи
!--------------ВОЗВРАЩАЕТ-------------------------
! ipar(1) - размерность
! ipar(2) - число мономов
! ipar(3) - число коэффициентов
! ipar(4) - число точек выдачи
! rpar(1) - требуемая относительная погрешность
! rpar(2) - требуемая абсолютная погрешность
! ln - линейность
integer,parameter :: fid=301,max_len=1000000
character(1),parameter::cc='%',em='^'
integer::eof,cn,ipar(:),i,ax
character(len=1000)::st,st1
character(*)::fpar(:),cfg
real::rpar(:),xx
logical :: ln
rpar=0.0
ipar=0
ln=.false.
fpar=em
! поиск ключ. слов в конф. файле
open(unit=fid,file=trim(cfg),status='old')
do
	read(fid,'(a)',iostat=eof) st 
	if(eof<0) then 
		exit
	else if(eof>0) then
		stop 'CONFIGURATION FILE ERROR'
	else
		continue
	end if
	cn=index(st,cc)
	select case(cn)
	case(0)
		st1=st 
	case(1)
		continue
	case default
		st1=st(1:cn-1)
	end select
	call tsmr_search(st1,rpar,fpar)
	st1=''
end do
close(fid)
! проверка погрешностей
if(rpar(1)==0.0 .and. rpar(2)/=0.0) then
	rpar(1)=rpar(2)
else if(rpar(1)/=0.0 .and. rpar(2)==0.0) then
	rpar(2)=rpar(1)
else if(rpar(1)==0.0 .and. rpar(2)==0.0) then
	stop 'ERROR IN REL_TOL AND ABS_TOL'
else
	continue
end if
! линейность
if(fpar(2)==em) ln=.true.
! размерность
open(unit=fid,file=trim(fpar(1)),status='old')
	read(fid,fmt=*,iostat=eof) (xx,i=1,max_len) 
	if(eof/=0) ipar(1)=i-1 
close(fid)
!число мономов
if(ln .eqv. .false.) then
	open(unit=fid,file=trim(fpar(2)),status='old')
		read(fid,fmt=*,iostat=eof) (xx,i=1,max_len) 
		if(eof/=0) ipar(2)=i-1
		if(mod(ipar(2),2)==0) then
            ipar(2)=ipar(2)/2
        else
            stop 'ERROR IN SCHEME'
        END IF
	close(fid)
end if
!число ненулевых коэффициентов
open(unit=fid,file=trim(fpar(3)),status='old')
	read(fid,fmt=*,iostat=eof) (xx,i=1,max_len) 
	if(eof/=0) ipar(3)=i-1
	if(mod(ipar(3),3)==0) then
            ipar(3)=ipar(3)/3
    else
        stop 'ERROR IN COEFFICIENTS'
    END IF
close(fid)
!число точек вывода
open(unit=fid,file=trim(fpar(6)),status='old')
	read(fid,fmt=*,iostat=eof) (xx,i=1,max_len)  
	if(eof/=0) ipar(4)=I-1
close(fid)
contains

subroutine tsmr_search(st,rpar,fpar)
! ПОИСК КЛЮЧЕВЫХ СЛОВ
! -------------ВХОДНЫЕ----------------------------
! st - строка конфигурационного файла
! rpar - массив с вещественными числами
! rpar - массив вещественных чисел
! fpar - массив строк
!--------------ВОЗВРАЩАЕТ-------------------------
! rpar(1) - требуемая относительная погрешность
! rpar(2) - требуемая абсолютная погрешность
! fpar(1) - файл с начальными данными
! fpar(2) - файл со схемой
! fpar(3) - файл с коэффициентами
! fpar(4) - файл с таблицей
! fpar(5) - файл для записи выдачи
! fpar(6) - файл с точками выдачи
integer,parameter::kn=9,skn=7
character(*),parameter::sc='='
character(len=1000)::st,key(9)
character(*)::fpar(:) 
integer::sp,i
real::rpar(:)
data key /'rel_tol','abs_tol',&
'init_cond','scheme','coeff',&
'table','out_file','integ_points', 'mass_planets'/  
sp=index(st,sc)+1
do i=1,kn
	if(index(st,trim(adjustl(key(i))))==0) then
		continue
	else
		select case(i)
		case(1)
			read(st(sp:),*) rpar(1)
		case(2)
			read(st(sp:),*) rpar(2)
		case default
			read(st(sp:),'(a)') fpar(i-(kn-skn))
		end select 
	end if
end do
end subroutine tsmr_search
end subroutine tsmr_cfg

end module tsmr