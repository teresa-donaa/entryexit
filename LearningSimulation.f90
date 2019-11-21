MODULE LearningSimulation
!
USE globals
USE QL_routines
USE omp_lib
USE ifport
!
! Computes Monte Carlo Q-Learning simulations
!
IMPLICIT NONE
!
CONTAINS
!
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!
    SUBROUTINE computeExperiment ( iExperiment, codExperiment, alpha, ExplorationParameters, delta )
    !
    ! Computes statistics for one model
    !
    IMPLICIT NONE
    !
    ! Declaring dummy variables
    !
    INTEGER, INTENT(IN) :: iExperiment, codExperiment
    REAL(8), DIMENSION(numAgents), INTENT(IN) :: alpha, delta
    REAL(8), DIMENSION(numExplorationParameters) :: ExplorationParameters
    !
    ! Declaring local variable
    !
    INTEGER :: idumIP, ivIP(32), iyIP, idum2IP, &
        idum, iv(32), iy, idum2, &
        idumQ, ivQ(32), iyQ, idum2Q, &
         idumEntry, ivEntry(32), iyEntry, idum2Entry
    INTEGER :: iIters, iItersFix, i, j, h, l, iSession, iItersInStrategy, convergedSession, numSessionsConverged
    INTEGER :: state, statePrime, stateFix, actionPrime
    INTEGER, DIMENSION(numStates,numAgents) :: strategy, strategyPrime, strategyFix
    INTEGER :: pPrime(numAgents), p(numAgents)
    INTEGER :: iAgent, iState, iPrice, jAgent
    INTEGER :: minIndexStrategies, maxIndexStrategies
    INTEGER(8) :: numSessions_I8
    REAL(8), DIMENSION(numStates,numPrices,numAgents) :: Q
    REAL(8) :: uIniPrice(numAgents,numSessions), uExploration(2,numAgents), uEntry
    REAL(8) :: u(2), eps(numAgents)
    REAL(8) :: newq, oldq
    REAL(8) :: meanTimeToConvergence, seTimeToConvergence, medianTimeToConvergence
    REAL(8) :: EntryProb, ExitProb
    CHARACTER(len = 25) :: QFileName
    CHARACTER(len = LengthFormatTotExperimentsPrint) :: iSessionsChar, codExperimentChar
    CHARACTER(len = 200) :: PTrajectoryFileName
    LOGICAL :: maskConverged(numSessions)
    !
    ! Beginning execution
    !
    ! Initializing various quantities
    !
    converged = 0    
    indexStrategies = 0
    indexLastState = 0
    timeToConvergence = 0.d0
    !
    EntryProb = DemandParameters(5+2*numAgents)
    ExitProb = DemandParameters(6+2*numAgents)
    !
    WRITE(codExperimentChar,'(I0.<LengthFormatTotExperimentsPrint>)') codExperiment
    !
    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Loop over numSessions
    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !
    ! Generating uIniPrice
    !
    idumIP = -1
    idum2IP = 123456789
    ivIP = 0
    iyIP = 0
    CALL generate_uIniPrice(uIniPrice,idumIP,ivIP,iyIP,idum2IP)  
    !
    ! Starting loop over sessions
    !
    !$ CALL OMP_SET_NUM_THREADS(numCores)
    !$omp parallel do &
    !$omp private(idum,iv,iy,idum2,idumQ,ivQ,iyQ,idum2Q,Q,maxValQ,idumEntry,ivEntry,iyEntry,idum2Entry, &
    !$omp   strategyPrime,pPrime,p,statePrime,actionPrime,iIters,iItersFix,iItersInStrategy,convergedSession, &
    !$omp   state,stateFix,strategy,strategyFix,eps,uExploration,u,uEntry,oldq,newq,iAgent,iState,iPrice,jAgent, &
    !$omp   QFileName,iSessionsChar) &
    !$omp firstprivate(numSessions,PI,delta,uIniPrice,ExplorationParameters,itersPerEpisode,alpha, &
    !$omp   EntryProb,ExitProb,itersInPerfMeasPeriod,maxIters,printQ,codExperimentChar)
    DO iSession = 1, numSessions
        !
        PRINT*, 'Session = ', iSession, ' started'
        !
        ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ! Learning phase
        ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !
        ! Initializing random number generators
        !
        idum = -iSession
        idum2 = 123456789
        iv = 0
        iy = 0
        !
        idumQ = -iSession
        idum2Q = 123456789
        ivQ = 0
        iyQ = 0
        !
        ! Initializing Q matrices
        !
        !$omp critical
        CALL initQMatrices(iSession,EntryProb,idumQ,ivQ,iyQ,idum2Q,PI,delta,Q,maxValQ,strategyPrime)
        !$omp end critical
        strategy = strategyPrime
        !
        ! Randomly initializing prices and state
        !
        CALL initState(uIniPrice(:,iSession),EntryProb,ExitProb,p,statePrime,actionPrime)
        state = statePrime
        !
        ! Loop
        !
        iIters = 0
        iItersInStrategy = 0
        convergedSession = -1
        eps = ExplorationParameters(:numAgents)
        IF (typeExplorationMechanism .EQ. 4) eps = ExplorationParameters(:numAgents)
        !
        DO 
            !
            ! Iterations counter
            !
            iIters = iIters+1
            !
            ! Generating exploration random numbers
            !            
            CALL generateUExploration(uExploration,idum,iv,iy,idum2)  
            CALL generateUEntry(uEntry,idumEntry,ivEntry,iyEntry,idum2Entry)  
            !
            ! Compute pPrime by balancing exploration vs. exploitation
            !
            CALL computePPrime(ExplorationParameters,uExploration,uEntry,EntryProb,ExitProb, &
                strategyPrime,state,iIters,pPrime,Q,eps)
            !
            ! Defining the new state
            !
            p = pPrime
            statePrime = computeStateNumber(p)
            actionPrime = computeActionNumber(pPrime)
            !
            ! Each agent collects his payoff and updates
            !
            DO iAgent = 1, numAgents
                !
                ! Q matrices and strategies update
                !
                oldq = Q(state,pPrime(iAgent),iAgent)
                newq = oldq+alpha(iAgent)*(PI(actionPrime,iAgent)+delta(iAgent)*maxValQ(statePrime,iAgent)-oldq)
                Q(state,pPrime(iAgent),iAgent) = newq
                IF (newq .GT. maxValQ(state,iAgent)) THEN
                    !
                    maxValQ(state,iAgent) = newq
                    IF (strategyPrime(state,iAgent) .NE. pPrime(iAgent)) strategyPrime(state,iAgent) = pPrime(iAgent)
                    !
                END IF
                IF ((newq .LT. maxValQ(state,iAgent)) .AND. (strategyPrime(state,iAgent) .EQ. pPrime(iAgent))) THEN
                    !
                    CALL MaxLocBreakTies(numPrices-1,Q(state,:numPrices-1,iAgent),idumQ,ivQ,iyQ,idum2Q, &
                        maxValQ(state,iAgent),strategyPrime(state,iAgent))
                    !
                END IF
                !
            END DO
!@SP
IF (ANY(strategyPrime(:,:numAgents-1) .EQ. numPrices)) THEN
    !
    PRINT*, 'Halt! Unadmissible price found!!!'
    PAUSE
    !
END IF
!@SP
            !
            ! Measuring performance
            !
            IF ((PerfMeasPeriodTime .LE. 0) .AND. ((iIters-1)/itersPerEpisode .GE. ABS(PerfMeasPeriodTime))) THEN
                !
                ! PerfMeasPeriodTime < 0: at convergence after mandatory training of PerfMeasPeriodTime years
                !
                IF (ALL(strategyPrime(state,:) .EQ. strategy(state,:))) THEN
                    !
                    iItersInStrategy = iItersInStrategy+1
                    !
                ELSE
                    !
                    iItersInStrategy = 1
                    !
                END IF
                !
            ELSE 
                !
                ! PerfMeasPeriodTime > 1: Starts right after the completion of year PerfMeasPeriodTime - 1
                !
                IF ((PerfMeasPeriodTime .GE. 1) .AND. ((iIters-1)/itersPerEpisode .GE. (PerfMeasPeriodTime-1))) THEN
                    !
                    iItersInStrategy = iItersInStrategy+1
                    !
                END IF
                !
            END IF
            !
            ! Check for convergence in strategy
            !
            IF (convergedSession .EQ. -1) THEN
                !
                ! Maximum number of iterations exceeded
                IF (iIters .GT. maxIters) THEN
                    !
                    convergedSession = 0
                    strategyFix = strategy
                    stateFix = state
                    iItersFix = iIters
                    !
                END IF
                !
                ! Prescribed number of iterations reached
                IF ((iItersInStrategy .EQ. itersInPerfMeasPeriod) .AND. (PerfMeasPeriodTime .GE. 1)) THEN
                    !
                    convergedSession = 0
                    strategyFix = strategy
                    stateFix = state
                    iItersFix = iIters
                    !
                END IF
                !
                ! Convergence in strategy reached
                IF ((iItersInStrategy .EQ. itersInPerfMeasPeriod) .AND. (PerfMeasPeriodTime .LE. 0)) THEN
                    !
                    convergedSession = 1
                    strategyFix = strategy
                    stateFix = state
                    iItersFix = iIters
                    !
                END IF
                !
            END IF
            !
            ! Check for loop exit criterion
            !
            IF (convergedSession .NE. -1) EXIT
            !
            ! If no convergence yet, update and iterate
            !
            strategy(state,:) = strategyPrime(state,:)
            state = statePrime
            !
            ! End of loop over iterations
            !
        END DO          
        !
        ! Write Q matrices to file
        !
        IF (printQ .EQ. 1) THEN
            !
            ! Open Q matrices output file
            !
            !$omp critical
            WRITE(iSessionsChar,'(I0.5)') iSession
            QFileName = 'Q_' // TRIM(codExperimentChar) // '_' // iSessionsChar // '.txt'
            !
            ! Write on Q matrices to file
            !
            OPEN(UNIT = iSession,FILE = QFileName,RECL = 10000)
            DO iAgent = 1, numAgents
                !
                DO iState = 1, numStates
                    !
                    WRITE(iSession,*) Q(iState,:,iAgent)
                    !
                END DO
                !
            END DO
            CLOSE(UNIT = iSession)
            !$omp end critical
            !
        END IF
        !
        ! Record results at convergence
        !
        converged(iSession) = convergedSession
        timeToConvergence(iSession) = DBLE(iItersFix-itersInPerfMeasPeriod)/itersPerEpisode
        indexLastState(:,iSession) = convertNumberBase(stateFix-1,numPrices,LengthStates)
        indexStrategies(:,iSession) = computeStrategyNumber(strategyFix)
        !
        IF (convergedSession .EQ. 1) PRINT*, 'Session = ', iSession, ' converged'
        IF (convergedSession .EQ. 0) PRINT*, 'Session = ', iSession, ' did not converge'
        !
        ! End of loop over sessions
        !
    END DO
    !$omp end parallel do
    !
    ! Print InfoExperiment file
    !
    OPEN(UNIT = 996,FILE = FileNameInfoExperiment,STATUS = "REPLACE")
    DO iSession = 1, numSessions
        !
        WRITE(996,*) iSession
        WRITE(996,*) converged(iSession)
        WRITE(996,*) timeToConvergence(iSession)
        WRITE(996,*) indexLastState(:,iSession)
        DO iState = 1, numStates
            !
            WRITE(996,*) (indexStrategies((iAgent-1)*numStates+iState,iSession), iAgent = 1, numAgents)
            !
        END DO
        !
    END DO
    CLOSE(UNIT = 996)
    !
    ! Prints the RES output file
    !
    numSessionsConverged = SUM(converged)
    maskConverged = (converged .EQ. 1)
    meanNashProfitIn = SUM(NashProfitsIn)/numAgents
    meanCoopProfitIn = SUM(CoopProfitsIn)/numAgents
    meanNashProfitOut = SUM(NashProfitsOut(1:numAgents-1))/DBLE(numAgents-1)
    meanCoopProfitOut = SUM(CoopProfitsOut(1:numAgents-1))/DBLE(numAgents-1)
    !
    ! Time to convergence
    !
    meanTimeToConvergence = SUM(timeToConvergence,MASK = maskConverged)/numSessionsConverged
    seTimeToConvergence = &
        SQRT(SUM(timeToConvergence**2,MASK = maskConverged)/numSessionsConverged-meanTimeToConvergence**2)
    numSessions_I8 = numSessions
    CALL SORTQQ(LOC(timeToConvergence),numSessions_I8,SRT$REAL8)
    medianTimeToConvergence = timeToConvergence(NINT(0.5d0*numSessions))
    !
    ! Print output
    !
    IF (iExperiment .EQ. 1) THEN
        !
        WRITE(10002,891) &
            (i, i = 1, numAgents), &
            (i, i = 1, numExplorationParameters), (i, i = 1, numAgents), &
            (i, (j, i, j = 1, numAgents), i = 1, numAgents), &
            (i, i = 1, numDemandParameters), &
            (i, i = 1, numAgents), (i, i = 1, numAgents), &
            (i, i = 1, numAgents), (i, i = 1, numAgents), &
            (i, i = 1, numAgents), (i, i = 1, numAgents), &
            (i, i = 1, numAgents), (i, i = 1, numAgents), &
            (i, i = 1, numAgents), (i, i = 1, numAgents), &
            (i, i = 1, numAgents), (i, i = 1, numAgents), &
            ((i, j, j = 1, numPrices), i = 1, numAgents)
891     FORMAT('Experiment ', &
            <numAgents>('    alpha', I1, ' '), &
            <numExplorationParameters>('     beta', I1, ' '), <numAgents>('    delta', I1, ' '), &
            <numAgents>('typeQini', I1, ' ', <numAgents>('par', I1, 'Qini', I1, ' ')), &
            <numDemandParameters>('  DemPar', I0.2, ' '), &
            <numAgents>(' NashPriceIn', I1, ' '), <numAgents>('NashPriceOut', I1, ' '), &
            <numAgents>(' CoopPriceIn', I1, ' '), <numAgents>('CoopPriceOut', I1, ' '), &
            <numAgents>(' NashProftIn', I1, ' '), <numAgents>('NashProftOut', I1, ' '), &
            <numAgents>(' CoopProftIn', I1, ' '), <numAgents>('CoopProftOut', I1, ' '), &
            <numAgents>(' NashMktShIn', I1, ' '), <numAgents>('NashMktShOut', I1, ' '), &
            <numAgents>(' CoopMktShIn', I1, ' '), <numAgents>('CoopMktShOut', I1, ' '), &
            <numAgents>(<numPrices>('Ag', I1, 'Price', I0.2, ' ')), &
            '   numConv     avgTTC      seTTC     medTTC ')
        !
    END IF
    !
    WRITE(10002,9911) codExperiment, &
        alpha, MExpl, delta, &
        (typeQInitialization(i), parQInitialization(i, :), i = 1, numAgents), &
        DemandParameters, &
        NashPricesIn, NashPricesOut, CoopPricesIn, CoopPricesOut, &
        NashProfitsIn, NashProfitsOut, CoopProfitsIn, CoopProfitsOut, &
        NashMarketSharesIn, NashMarketSharesOut, CoopMarketSharesIn, CoopMarketSharesOut, &
        (PricesGrids(:,i), i = 1, numAgents), &
        numSessionsConverged, meanTimeToConvergence, seTimeToConvergence, medianTimeToConvergence
9911 FORMAT(I5, 1X, &
        <numAgents>(F10.5, 1X), <numExplorationParameters>(F10.5, 1X), <numAgents>(F10.5, 1X), &
        <numAgents>(A9, 1X, <numAgents>(F9.2, 1X)), &
        <numDemandParameters>(F10.5, 1X), &
        <6*numAgents+6*numAgents>(F13.3, 1X), &
        <numPrices*numAgents>(F10.3, 1X), &
        I10, 1X, <3>(F10.2, 1X))
    !
    ! Ending execution and returning control
    !
    END SUBROUTINE computeExperiment
!
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!
    SUBROUTINE computePPrime ( ExplorationParameters, uExploration, uEntry, EntryProb, ExitProb, &
        strategyPrime, state, iIters, pPrime, Q, eps )
    !
    ! Computes pPrime by balancing exploration vs. exploitation
    !
    IMPLICIT NONE
    !
    ! Declaring dummy variables
    !
    REAL(8), INTENT(IN) :: ExplorationParameters(numExplorationParameters)
    REAL(8), INTENT(IN) :: uExploration(2,numAgents)
    REAL(8), INTENT(IN) :: uEntry
    REAL(8), INTENT(IN) :: EntryProb, ExitProb
    INTEGER, INTENT(IN) :: strategyPrime(numStates,numAgents)
    INTEGER, INTENT(IN) :: state, iIters
    INTEGER, INTENT(OUT) :: pPrime(numAgents)
    REAL(8), INTENT(IN) :: Q(numStates,numPrices,numAgents)
    REAL(8), INTENT(INOUT) :: eps(numAgents)
    !
    ! Declaring local variables
    !
    INTEGER :: iAgent, iPrice
    REAL(8) :: u(2), maxQ
    REAL(8) :: probs(numPrices)
    !
    ! Beginning execution
    !
    ! 1. Greedy with probability 1-epsilon, with constant epsilon
    !
    IF (typeExplorationMechanism .EQ. 1) THEN
        !
        DO iAgent = 1, numAgents-1
            !
            eps(iAgent) = ExplorationParameters(iAgent)
            u = uExploration(:,iAgent)
            IF (u(1) .LE. eps(iAgent)) THEN
                !
                pPrime(iAgent) = 1+INT((numPrices-1)*u(2))
                !
            ELSE
                !
                pPrime(iAgent) = strategyPrime(state,iAgent)
                !
            END IF
            !
        END DO
        u = uExploration(:,numAgents)
        IF (MOD(state,numPrices) .EQ. 0) THEN       ! Last agent is currently OUT
            !
            IF (uEntry .LE. EntryProb) THEN           
                !
                IF (u(1) .LE. eps(numAgents)) THEN
                    !
                    pPrime(numAgents) = 1+INT((numPrices-1)*u(2))
                    !
                ELSE
                    !
                    pPrime(numAgents) = strategyPrime(state,numAgents)
                    !
                END IF
                !
            ELSE 
                ! 
                pPrime(numAgents) = numPrices
                !
            END IF
            !
        ELSE                                        ! Last agent is currently IN
            !
            IF (uEntry .LE. ExitProb) THEN           
                !
                pPrime(numAgents) = numPrices
                !
            ELSE 
                ! 
                IF (u(1) .LE. eps(numAgents)) THEN
                    !
                    pPrime(numAgents) = 1+INT((numPrices-1)*u(2))
                    !
                ELSE
                    !
                    pPrime(numAgents) = strategyPrime(state,numAgents)
                    !
                END IF
                !
            END IF
            !
        END IF
        !
    END IF
    !
    ! 2. & 3. Greedy with probability 1-epsilon, with exponentially decreasing epsilon
    !
    IF ((typeExplorationMechanism .EQ. 2) .OR. (typeExplorationMechanism .EQ. 3)) THEN
        !
        DO iAgent = 1, numAgents-1
            !
            IF (MExpl(numAgents+iAgent) .LT. 0.d0) THEN
                !
                pPrime(iAgent) = strategyPrime(state,iAgent)
                !
            ELSE
                !
                u = uExploration(:,iAgent)
                IF (u(1) .LE. eps(iAgent)) THEN
                    !
                    pPrime(iAgent) = 1+INT((numPrices-1)*u(2))
                    !
                ELSE
                    !
                    pPrime(iAgent) = strategyPrime(state,iAgent)
                    !
                END IF
                eps(iAgent) = eps(iAgent)*ExplorationParameters(numAgents+iAgent)
                !
            END IF
            !
        END DO
        IF (MExpl(2*numAgents) .LT. 0.d0) THEN
            !
            pPrime(numAgents) = strategyPrime(state,numAgents)
            !
        ELSE
            !
            u = uExploration(:,numAgents)
            IF (MOD(state,numPrices) .EQ. 0) THEN       ! Last agent is currently OUT
                !
                IF (uEntry .LE. EntryProb) THEN           
                    !
                    IF (u(1) .LE. eps(numAgents)) THEN
                        !
                        pPrime(numAgents) = 1+INT((numPrices-1)*u(2))
                        !
                    ELSE
                        !
                        pPrime(numAgents) = strategyPrime(state,numAgents)
                        !
                    END IF
                    !
                ELSE 
                    ! 
                    pPrime(numAgents) = numPrices
                    !
                END IF
                !
            ELSE                                        ! Last agent is currently IN
                !
                IF (uEntry .LE. ExitProb) THEN           
                    !
                    pPrime(numAgents) = numPrices
                    !
                ELSE 
                    ! 
                    IF (u(1) .LE. eps(numAgents)) THEN
                        !
                        pPrime(numAgents) = 1+INT((numPrices-1)*u(2))
                        !
                    ELSE
                        !
                        pPrime(numAgents) = strategyPrime(state,numAgents)
                        !
                    END IF
                    !
                END IF
                !
            END IF
            eps(numAgents) = eps(numAgents)*ExplorationParameters(2*numAgents)
            !
        END IF
        !
    END IF
    !
    ! 4. Boltzmann exploration
    !
    IF (typeExplorationMechanism .EQ. 4) THEN
        !
        DO iAgent = 1, numAgents-1
            !
            maxQ = MAXVAL(Q(state,:numPrices-1,iAgent))
            probs(:numPrices-1) = EXP((Q(state,:numPrices-1,iAgent)-maxQ)/eps(iAgent))
            u = uExploration(:,iAgent)*SUM(probs(:numPrices-1))
            DO iPrice = 1, numPrices-1
                !
                IF (iPrice .GT. 1) probs(iPrice) = probs(iPrice)+probs(iPrice-1)
                IF (u(1) .LE. probs(iPrice)) THEN
                    !
                    pPrime(iAgent) = iPrice
                    EXIT
                    !
                END IF
                !
            END DO
            eps(iAgent) = eps(iAgent)*ExplorationParameters(numAgents+iAgent)
            !
        END DO
        IF (MOD(state,numPrices) .EQ. 0) THEN       ! Last agent is currently OUT
            !
            IF (uEntry .LE. EntryProb) THEN           
                !
                maxQ = MAXVAL(Q(state,:numPrices-1,numAgents))
                probs(:numPrices-1) = EXP((Q(state,:numPrices-1,numAgents)-maxQ)/eps(numAgents))
                u = uExploration(:,numAgents)*SUM(probs(:numPrices-1))
                DO iPrice = 1, numPrices-1
                    !
                    IF (iPrice .GT. 1) probs(iPrice) = probs(iPrice)+probs(iPrice-1)
                    IF (u(1) .LE. probs(iPrice)) THEN
                        !
                        pPrime(numAgents) = iPrice
                        EXIT
                        !
                    END IF
                    !
                END DO
                !
            ELSE 
                ! 
                pPrime(numAgents) = numPrices
                !
            END IF
            !
        ELSE                                        ! Last agent is currently IN
            !
            IF (uEntry .LE. ExitProb) THEN           
                !
                pPrime(numAgents) = numPrices
                !
            ELSE 
                ! 
                maxQ = MAXVAL(Q(state,:numPrices-1,numAgents))
                probs(:numPrices-1) = EXP((Q(state,:numPrices-1,numAgents)-maxQ)/eps(numAgents))
                u = uExploration(:,numAgents)*SUM(probs(:numPrices-1))
                DO iPrice = 1, numPrices-1
                    !
                    IF (iPrice .GT. 1) probs(iPrice) = probs(iPrice)+probs(iPrice-1)
                    IF (u(1) .LE. probs(iPrice)) THEN
                        !
                        pPrime(numAgents) = iPrice
                        EXIT
                        !
                    END IF
                    !
                END DO
                !
            END IF
            !
        END IF
        eps(numAgents) = eps(numAgents)*ExplorationParameters(2*numAgents)
        !
    END IF
    !
    ! Ending execution and returning control
    !
    END SUBROUTINE computePPrime
!
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!
END MODULE LearningSimulation
