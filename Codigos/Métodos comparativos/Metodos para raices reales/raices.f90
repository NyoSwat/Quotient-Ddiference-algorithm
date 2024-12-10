
!PROGRAMA UTILIZADO PARA CALCULAR LAS RAICES POR LOS MÉTODOS VISTOS EN CLASE
PROGRAM Raices
IMPLICIT NONE
REAL(8) :: a,b
REAL(8),PARAMETER :: tolerancia=0.0001
INTEGER,PARAMETER :: maxiter=100
!En el parametro metodo se selecciona lo que se desea realizar:      
!1: Graficar la funcion.                                            
!2: Metodo Biseccion.                                               
!3: Metodo Punto fijo sistematico.                                  
!4: Metodo Newton-Raphson.
INTEGER,PARAMETER :: metodo=2                                   

OPEN(UNIT=1,FILE='resolucion.dat',STATUS='REPLACE')
a=1. !Extremo inferior
b=3. !Extremo superior
CALL verificaextremos(a,b)

SELECT CASE(metodo) 
    CASE(1)
        CALL grafica(a,b)
    CASE(2)
        CALL biseccion(a,b)
        WRITE(*,*)'La resolucion se encuentra en archivo resolucion.dat'
    CASE(3)
        CALL puntofijo(a,b)
        WRITE(*,*)'La resolucion se encuentra en archivo resolucion.dat'
    CASE(4)
        CALL newtonraphson(a)
        WRITE(*,*)'La resolucion se encuentra en archivo resolucion.dat'
END SELECT

CLOSE(UNIT=1,STATUS='KEEP')

CONTAINS

SUBROUTINE verificaextremos(a,b)!VERIFICA QUE EL EXTREMO INFERIOR SEA MENOR AL SUPERIOR

REAL(8) :: a,b

IF (a>b) THEN
    WRITE (*,*)'a>b! Se ingresaron incorrectamente los extremos.'
END IF

END SUBROUTINE verificaextremos

SUBROUTINE grafica(a,b)

REAL(8) :: a,b

OPEN(UNIT=2,FILE='grafica.p',STATUS='REPLACE')

 WRITE(2,*)'set xrange[',a,':',b,' ]'
 WRITE(2,*)'unset log                              # quita la escala logaritmica (si la hubiera)'
 WRITE(2,*)'unset label                            # quita los titulos anteriores'
 WRITE(2,*)'set xtic auto                          # establece automaticamente las divisiones del eje x'
 WRITE(2,*)'set ytic auto                          # establece automaticamente las divisiones del eje y'
 WRITE(2,*)'set grid'
 WRITE(2,*)'set title " Grafico de f(x)=0 para hallar raiz "'
 WRITE(2,*)'set xlabel "x"'
 WRITE(2,*)'set ylabel "y"'
 WRITE(2,'(a)',ADVANCE='NO')'plot '
 !En el siguiente renglon se escribe la funcion a graficar.
 WRITE(2,'(a)',ADVANCE='NO')'x*exp(x)-1'
 WRITE(2,'(a)',ADVANCE='NO')' with lines'
 CALL SYSTEM (" C:\gnuplot\bin\gnuplot.exe -persist grafica.p")

CLOSE(UNIT=2,STATUS='KEEP')

END SUBROUTINE grafica

FUNCTION f(x) !Funcion
REAL(8) x,f
       f=x**3 -4.*x**2 +x+6.
END FUNCTION f

FUNCTION fprima(x) !Derivada de la funcion
REAL(8) x,fprima

   fprima= 3.*x**2 - 8.*x +1

END FUNCTION fprima

FUNCTION maxderivada(a,b) !Maximo valor de la derivada

REAL(8) :: q,a,b,paso,maxderivada,t

q=a
paso=0.001
maxderivada=0.
DO WHILE (q<=b)
    t=fprima(q)
    IF (abs(t)>abs(maxderivada)) THEN
        maxderivada=t
    END IF
    q=q+paso
END DO
END FUNCTION maxderivada

SUBROUTINE biseccion(a,b)

INTEGER :: i,n
REAL(8) :: a,b,a2,b2,m,errorenx,erroreny,raiz, xant, errorpaso

WRITE(1,'(A/)')' RESOLUCION MEDIANTE EL METODO DE BISECCION '
WRITE(1,*)
n=FLOOR((LOG(ABS(b-a)/tolerancia))/(LOG(2.0))+0.5)
WRITE(*,'(A,I4,A/)')'Se realizaran aproximadamente ',n,' iteraciones.'
a2=a
b2=b
erroreny=99999.0
i=0
m=(a2+b2)/2.0
errorpaso=1000.0
DO WHILE (tolerancia<ABS(erroreny))
!   DO WHILE (tolerancia < ABS(errorpaso))
    i=i+1
    xant=m
    errorenx=(abs(b2-a2)/2.0)
    erroreny=f((a2+b2)/2.0)
    WRITE(1,*)'------------------------------------------------------------------'
    WRITE(1,*)'Iteración: ',i
    WRITE(1,'(A,F20.7,A,F20.7,A)')'Se analiza intervalo: [ ',a2,' , ',b2,' ]'
    WRITE(1,'(A,F20.7)')'El error en x es: ',errorenx
    WRITE(1,'(A,F20.7)')'El error en y es: ',erroreny
    WRITE(1,*)
    IF (f(a2)*f(m)<0) THEN
        b2=m
        WRITE(1,*)'B pasara a ser: ',b2
        WRITE(1,*)'------------------------------------------------------------------'
    ELSE
        a2=m
        WRITE(1,*)'A pasara a ser: ',a2
        WRITE(1,*)'------------------------------------------------------------------'
    END IF
    m=(a2+b2)/2.0
    errorpaso= xant-m
END DO
IF ((ABS(erroreny)<tolerancia)) THEN
    errorenx=(abs(b2-a2)/2.0)
    erroreny=f(m)
    raiz=m
    WRITE(1,*) '------------------------------------------------------------------'
    WRITE(1,'(A,F20.7)')'La raiz encontrada es: ',raiz
    WRITE(1,'(A,F20.7)')'El error en x es: ',errorenx
    WRITE(1,'(A,F20.7)')'El error en y es: ',erroreny
ELSE
    WRITE(1,'(A/)')'No hay raices reales en el intervalo o se trata de una SINGULARIDAD'
END IF
WRITE(*,'(A,F20.7/)')'La raiz encontrada es: ',raiz

END SUBROUTINE biseccion

SUBROUTINE puntofijo(a,b)

REAL(8) :: a,b,x,lamda,error, xant, errorenx
INTEGER :: iter

WRITE(1,'(A/)')' RESOLUCION MEDIANTE PUNTO FIJO SISTEMATICO '
error=(tolerancia*2.0)
errorenx=(tolerancia*2.0)
iter=1
x=a
lamda=(1.0/maxderivada(a,b))
DO WHILE (error>tolerancia) 
!  DO WHILE (errorenx>tolerancia) 
    WRITE(1,*)'---------------------------------------'
    xant=x
    x=x-(lamda*f(x))
    WRITE(1,'(A,I4)')'Iteracion: ',iter
    WRITE(1,'(A,F20.7)')'Valor de x: ',x
    error=abs(f(x))
    WRITE(1,'(A,F20.7)')'Error o distancia a y: ',error
    iter=iter+1
    errorenx= abs(xant-x)
    WRITE(1,'(A,F20.7)')'Error en x: ',errorenx
END DO
WRITE(1,*)'---------------------------------------'
WRITE(1,'(A,F20.7)')'La raiz encontrada es: ',x
WRITE(*,'(A,F20.7)')'La raiz encontrada es: ',x

END SUBROUTINE puntofijo

SUBROUTINE newtonraphson (a)

REAL(8) :: a,x,error, errorenx, xant
INTEGER :: iter

WRITE(1,'(A/)') ' RESOLUCION MEDIANTE EL METODO DE NEWTONRAPSHON '
error=(tolerancia*2.0)
errorenx=(tolerancia*2.0)
iter=1
WRITE(*,*)'Ingresar el valor de "x0".'
WRITE(*,*)'En caso de no tenerlo, puede colocar el extremo "a".'
WRITE(*,*)'Recordar que el extremos "a" vale: ',a
READ(*,*)x
DO WHILE (error>=tolerancia)  
!   DO WHILE (errorenx> tolerancia)
    WRITE(1,*)'---------------------------------------'
    WRITE(1,*)
    xant=x
    x=x-(f(x)/fprima(x))
    WRITE(1,'(A,I4)')'Iteracion: ',iter
    WRITE(1,'(A,F20.7)')'Valor de x: ',x
    error=abs(f(x))
    WRITE(1,'(A,F20.7)')'Error o distancia a y: ',error
    iter=iter+1
    errorenx= abs(x-xant)
    WRITE(1,'(A,F20.7)')'Error en x : ',errorenx
END DO
WRITE(1,*)'---------------------------------------'
WRITE(1,*)
WRITE(1,'(A,F20.7/)')'La raiz encontrada es: ',x
WRITE(*,'(A,F20.7/)')'La raiz encontrada es: ',x

END SUBROUTINE newtonraphson

END PROGRAM
