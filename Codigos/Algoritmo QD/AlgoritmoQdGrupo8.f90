PROGRAM ALGORITMO_QD

    IMPLICIT NONE 
    INTEGER N, i, iteracionMax,opcion
    REAL(8), ALLOCATABLE, DIMENSION(:) :: P, PAux1, PAux2
    COMPLEX(8), ALLOCATABLE, DIMENSION(:) :: Raices
    REAL(8) tolerancia
    REAL(8) c
  
    OPEN (1,FILE ="coeficientes.txt")
    READ(1,*) N
    ALLOCATE(P(0:N), PAux1(0:N), Raices(1:N),PAux2(0:N))
    do i = 0, N
        READ(1,*) P(i)
    end do
    
        
    PAux1 = P   !Polinomio auxiliar en caso de usar la traslacion. 
    WRITE(*,*)
    WRITE(*,*)'           Quotient-Difference -algorithm ' 
    WRITE(*,*) 
    ! WRITE(*,*) 'El polinomio original es el siguiente:'
    WRITE(*,*) 'El polinomio ingresado es el siguiente:' 
    WRITE(*,*) 
    ! WRITE(*,'(A)', ADVANCE='NO') ' P(x) = ' (dentro de la subrutina)
    CALL Muestra_Polinomio(P, N)
    WRITE(*,*)
    WRITE(*,*) 'Realizar una transformacion sobre la indeterminada al polinomio? 1 = SI, 2 = NO'
    READ(*,*) opcion

    IF (opcion == 1) THEN
        WRITE(*,*) 'Ingrese el factor c para efectuar la traslacion.'
            READ(*,*) c
        DO WHILE(c == 0)
            WRITE(*,*) 'Factor invalido.'
            WRITE(*,*) 'Ingrese el factor c nuevamente.'
            READ(*,*) c     
            END DO
        CALL Translacion_x(PAux1, N, c)  
        PAux2 = PAux1
    ELSE
        PAux2 = P
        
    END IF
                                                        
    tolerancia = 0.0001
    iteracionMax = 3     !Carga de la toleranciaerancia y la cantidad maxima de iteraciones. 

    CALL QD(PAux2, Raices, N, tolerancia, iteracionMax)
    CALL grafica_polinomio(Raices,N,P)
    CALL Muestra_Raices(Raices, N, tolerancia)
    DEALLOCATE(P, PAux1,PAux2, Raices)

CONTAINS
!----------------------------------------------------------------------------------------------------------------

SUBROUTINE Muestra_Polinomio(P, N) 
    REAL(8), DIMENSION(0:N),INTENT(INOUT) :: P
    INTEGER,INTENT(IN) :: N
    INTEGER i
        WRITE(*,'(A)', ADVANCE='NO') ' P(x) = ' 
        DO i=N, 1, -1
            IF (P(i) /= 0 ) THEN
                WRITE(*,'(F10.3, A4, I1, A3)', ADVANCE='NO') P(i), ' x**', i, ' + ' !Escribimos en pantalla todos los términos acompañados de x
            END IF
        END DO
        WRITE(*,'(F10.3, A, I1, A)', ADVANCE='NO') P(i) !Escribimos en pantalla el término independiente
        WRITE(*,*)

END SUBROUTINE Muestra_Polinomio

!----------------------------------------------------------------------------------------------------

SUBROUTINE QD(P, Raices, N, tolerancia, iteracionMax)

    REAL(8), DIMENSION(0:N), INTENT(INOUT) :: P
    REAL(8), DIMENSION(0:N) :: e
    REAL(8), DIMENSION(1:N) :: q , q_ant
    REAL(8), INTENT(IN) :: tolerancia
    REAL(8) u, v, Error
    COMPLEX(8), DIMENSION(1:N), INTENT(INOUT) :: Raices
    COMPLEX(8) c1, c2
    INTEGER, INTENT(IN) :: N, iteracionMax
    INTEGER i, iter

    q = 0.
    q(1) = -P(N-1) / P(N)
    e(0) = 0.
    e(N) = 0.
    
    DO i=1, N-1
        e(i) = P(N-i-1) / P(N-i)
    END DO
    
    iter = 0
    Error = 2. * tolerancia
    
    OPEN(2, FILE='Tabla.txt')
    WRITE(*,*)
    WRITE(2,'(A)', ADVANCE='NO') 'Iteración        e0         '   
    
    DO i=1, N 
        WRITE(2,'(A, I1, A, I1, A)', ADVANCE='NO') 'q', i,'         e', i,'         '
    END DO !Graba en un txt la tabla con los valores de q y e en la iteracion 0.
    WRITE(2,*)    

    DO WHILE (error >= tolerancia .AND. iter <= iteracionMax) 
        WRITE(*,'(I5, A)', ADVANCE='NO') iter, '      '
        DO i=1, N
            WRITE(*,'(F22.5)', ADVANCE='NO')  q(i)
        END DO
        WRITE(*,*) !Muestra por pantalla los valores de q y e en cada iteracion distinta de 0.
        DO i=0, N
            WRITE(*,'(F22.5)', ADVANCE='NO') e(i)
        END DO
        WRITE(*,*)
        WRITE(2,'(I5, A)', ADVANCE='NO') iter, '      '
        DO i=1, N
            WRITE(2,'(F22.5)', ADVANCE='NO')  q(i)
        END DO
        WRITE(2,*) !Graba en un txt los valores de q y e en cada iteracion distinta de 0.
        DO i=0, N
            WRITE(2,'(F22.5)', ADVANCE='NO') e(i)
        END DO
        WRITE(2,*)
        
        q_ant = q
        DO i=1, N
            q(i) = e(i) - e(i-1) + q(i)
        END DO !Genero los nuevos valores de q y e. 
        DO i=1, N-1
            e(i) = (q(i+1) / q(i)) * e(i)
        END DO
        Error = maxval(abs(e))
        iter = iter+1        
    END DO
    WRITE(*,*)
    DO i=0, N-1
        IF ((abs(e(i)) > tolerancia)) THEN !Si el error fluctua
            u = q(i) + q(i+1)
            v = -(q_ant(i) * q(i+1))
            CALL Bairstow (P, N, u, v, tolerancia) !Calculamos, de ser necesario las raices complejas o co-modulares.
            CALL Resolvente(-u, -v, c1, c2)
            Raices(i) = c1 !Guardamos todas las raices, ya sean reales o complejas en un unico vector
            Raices(i+1) = c2
        ELSE
            Raices(i+1) = DCMPLX(q(i+1), 0.)
        END IF
    END DO
    CLOSE(2)

END SUBROUTINE QD

!------------------------------------------------------------------------------------------------

SUBROUTINE Bairstow(A, N, u, v, tolerancia)

    REAL(8), DIMENSION(0:N), INTENT(INOUT) :: A
    REAL(8), INTENT(INOUT) :: u, v
    REAL(8), INTENT(IN) :: tolerancia
    REAL(8) q, q_ant1, q_ant2 , p, p_ant1, p_ant2, p_ant3, h, k, Error
    INTEGER, INTENT(IN) :: N
    INTEGER i

    Error = 2. * tolerancia
    DO WHILE(Error >= tolerancia)
        
        q_ant1 = 0.
        q_ant2 = 0.
        p_ant1 = 0. !valores inciales para la aplicacion de Bairstow
        p_ant2 = 0.
        DO i=N, 1, -1
            q = A(i) + u * q_ant1 + v * q_ant2
            q_ant2 = q_ant1
            q_ant1 = q
            p = q + u * p_ant1 + v * p_ant2   !Calculo de los valores de q y p para i=1... N
            p_ant3 = p_ant2
            p_ant2 = p_ant1
            p_ant1 = p
        END DO
        q = A(0) + u * q_ant1 + v * q_ant2  !Calculo del valor de q para i=0
        
        h = (q * p_ant3 - q_ant1 * p_ant2) / (p_ant2**2. - p_ant1 * p_ant3)
        k = (q_ant1 * p_ant1 - q * p_ant2) / (p_ant2**2. - p_ant1 * p_ant3)
        
        u = u + h
        v = v + k
        IF(abs(q) > abs(q_ant1))THEN
            Error = abs(q)
        ELSE
            Error= abs(q_ant1)
        END IF
    END DO
        
END SUBROUTINE Bairstow

!-------------------------------------------------------------------------------------------------

SUBROUTINE Resolvente(b, c, c1, c2)

    REAL(8) a
    REAL(8), INTENT(IN) :: b, c
    COMPLEX(8) disc 
    COMPLEX(8), INTENT(OUT) :: c1, c2

    a = 1.
    disc = b**2 - 4. * a * c
    c1 = (-b + sqrt(disc)) / (2. * a)
    c2 = (-b - sqrt(disc)) / (2. * a)
    WRITE(*,*)'Raices calculadas por el metodo de Bairstow'
    WRITE(*,*)c1,c2
    WRITE(*,*)

END SUBROUTINE Resolvente

!--------------------------------------------------------------------------------------------

SUBROUTINE Muestra_Raices(Raices, N, tolerancia)

 !REAL(8), DIMENSION(0:N),INTENT(INOUT) :: P  
 COMPLEX(8), DIMENSION(1:N),INTENT(INOUT) :: Raices
 REAL(8),INTENT(IN) :: tolerancia
 INTEGER,INTENT(INOUT) :: N
 INTEGER i

    WRITE(*,'(A)', ADVANCE='NO') 'Raices del Polinomio P(x) = '
    !DO i=N, 1, -1
    !    WRITE(*,'(F10.5, A, I2, A)', ADVANCE='NO') P(i), '* x**', i, " +" 
    !END DO
    !WRITE(*,'(F10.5)') P(0)                                                    
    i=1
    DO WHILE (i <= N)
        IF (IMAG(Raices(i)) == 0.) THEN
            WRITE(*,'(A, I2, A, F10.5)') ' x', i, ' = ', REAL(Raices(i))
        ELSE IF (abs(REAL(Raices(i))) <= tolerancia ) THEN
            WRITE(*,'(A, I2, A, F10.5, A, F10.5, A)') ' x', i, ' = ', IMAG(Raices(i)), ' i'
            i = i+1
            WRITE(*,'(A, I2, A, F10.5, A, F10.5, A)') ' x', i, ' = ', IMAG(Raices(i)), ' i'  !Muestro por pantalla las raices con un formato correspondiente
        ELSE 
            WRITE(*,'(A, I2, A, F10.5, A, F10.5, A)') ' x', i, ' = ', REAL(Raices(i)), ' + ', IMAG(Raices(i)), ' i' 
            i = i+1
            WRITE(*,'(A, I2, A, F10.5, A, F10.5, A)') ' x', i, ' = ', REAL(Raices(i)), ' - ', abs(IMAG(Raices(i))), ' i' 
        END IF
        i = i+1
    END DO

END SUBROUTINE Muestra_Raices

!-------------------------------------------------------------------------------------------------------------------

SUBROUTINE Translacion_x(P, N, c)

    REAL(8), DIMENSION(0:N), INTENT(INOUT) :: P
    REAL(8), ALLOCATABLE, DIMENSION(:) :: PAux1
    REAL(8), INTENT(IN) :: c
    INTEGER, INTENT(IN) :: N
    INTEGER i, tope

    ALLOCATE(PAux1(0:N))
    PAux1 = P
    tope = N
    
    DO i=0, N
        P(i) = Evalua_Polinomio(PAux1, tope, c) / Factorial(i)
        CALL Derivada(PAux1, tope)
    END DO
    DEALLOCATE(PAux1)

END SUBROUTINE Translacion_x

!-----------------------------------------------------------------------

SUBROUTINE Derivada(P, N)

    REAL(8), DIMENSION(0:N), INTENT(INOUT) :: P
    INTEGER, INTENT(INOUT) :: N
    INTEGER i

    DO i=0, N-1
        P(i) = P(i+1) * (i+1)
    END DO
    N = N - 1

END SUBROUTINE Derivada

!-----------------------------------------------------------------------

FUNCTION Evalua_Polinomio(P, N, x)

    REAL(8), DIMENSION(0:N), INTENT(IN) :: P
    REAL(8), INTENT(IN) :: x 
    REAL(8) Suma, Evalua_Polinomio
    INTEGER, INTENT(IN) :: N
    INTEGER i

    Suma = 0.
    DO i=N, 1, -1
        Suma = (Suma + P(i)) * x
    END DO
    Evalua_Polinomio = Suma + P(0)

END FUNCTION Evalua_Polinomio

!-----------------------------------------------------------------------

FUNCTION Factorial(x)
    REAL(8) Factorial
    INTEGER, INTENT(IN) :: x
    INTEGER i, aux

    IF (x == 0 .OR. x==1) THEN
        Factorial = 1.
    ELSE 
        aux = 1
        DO i=2, x
           aux = aux * i
        END DO
    Factorial = aux
    END IF

END FUNCTION Factorial

!------------------------------------------------------------------------ 

FUNCTION f(x,P,N)

    REAL(8) x,f
    INTEGER, INTENT(IN) :: N
    REAL(8), DIMENSION(0:N), INTENT(in):: P
    do i = 0, N
        f = f + P(i)*(x**i)
    end do

END FUNCTION

!----------------------------------------------------------------------------------------------------------------

SUBROUTINE Grafica_Polinomio(Raices,N,P)
    
    REAL (8) x
    REAL(8), PARAMETER :: paso_grafico = 0.05
    COMPLEX(8), DIMENSION(1:N),INTENT(INOUT) :: Raices
    INTEGER, INTENT(IN) :: N
    REAL(8), DIMENSION(0:N), INTENT(IN) :: P
    INTEGER i
    
    OPEN (UNIT=2, FILE ='RaicesReales.dat', STATUS='REPLACE')
    OPEN (UNIT=5, FILE ='RaicesComplejas.dat', STATUS='REPLACE')

    DO i=1,N
        IF (imag(Raices(i)) .ne. 0) THEN
            WRITE (5, '(2F22.5)')  Raices(i),0.0
        ELSE                                       !Solo se grafican las raices reales
            WRITE (2, '(2F22.5)')  Raices(i),0.0
        END IF 
    END DO
        
    CLOSE(2)
    
    OPEN (UNIT=3, FILE ='Polinomio.dat', STATUS='REPLACE') 
    x = -10  !Limite inferior del grafico de la funcion
    DO WHILE (x <= 10)  !Limite superior del grafico de la funcion
        WRITE (3, '(2F22.5)')  x , f(x,P,N)   
        x = x + paso_grafico
    END DO
    CLOSE(3)

    !script para gnuplot
    CALL SYSTEM('gnuplot -persist script.p') !reemplazar 'gnuplot' por el camino al ejecutable en caso de no funcionar correctamente
    END SUBROUTINE Grafica_Polinomio
    
!-----------------------------------------------------------------------

END PROGRAM ALGORITMO_QD

