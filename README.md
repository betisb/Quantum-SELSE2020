# Quantum 9X4-qubit error correction algorithm (Q9X4-ECC) 

###### this tutorial intent to practice Quantum error correction with additional extra twist in Q# code over Shor’s 9-qubit error correcting code (ECC) algorithm Shor [1995]. 

* To begin, first prepare this notebook for execution (if you skip this step, you'll get "Syntax does not match any known patterns" error when you try to execute Q# code in the next cells; if you skip the second step, you'll get "Invalid test name" error):

        In [ ]: %package Microsoft.Quantum.Katas::0.10.1911.1607

The package versions in the output of the cell above should always match. If you are running the Notebooks locally and the versions do not match, please install the IQ# version that matches the version of the Microsoft.Quantum.Katas package.

    In [ ]: %workspace reload
.




    function ToString_Bitmask (bits : Int) : String {
        let b1 = bits / 4;
        let b2 = (bits / 2) % 2;
        let b3 = bits % 2;
        let b4 = bits % 3;
        let b5 = bits % 5;
        let b6 = bits % 7;
        let b7 = bits % 9;
        let b8 = bits % 11;
        let b9 = bits % 13;
        return $"{b1}{b2}{b3}{b4}{b5}{b6}{b7}{b8}{b9}";
    }
    operation StatePrep_Bitmask (qs : Qubit[], bits : Int) : Unit
    is Adj {
        
        if (bits / 4 == 1) {
            X(qs[0]);
        }
            
        if ((bits / 2) % 2 == 1) {
            X(qs[1]);
        }
            
        if (bits % 2 == 1) {
            X(qs[2]);
        }
        if (bits % 3 == 1){
            X(qs[3]);
        }
        if (bits % 5 == 1){
            X(qs[4]);
        }
        if (bits % 7 == 1){
            X(qs[5]);
        }
        if (bits % 9 == 1){
            X(qs[6]);
        }
        if (bits % 11 == 1){
            X(qs[7]);
        }
        if (bits % 13 == 1){
            X(qs[8]);
        }
    }
        
    // same Function as Katas Kit
    
    function FindFirstDiff_Reference (bits1 : Int[], bits2 : Int[]) : Int {
            mutable firstDiff = -1;
            for (i in 0 .. Length(bits1) - 1) {
                if (bits1[i] != bits2[i] and firstDiff == -1) {
                    set firstDiff = i;
                }
            }
            return firstDiff;
        }
        
    operation StatePrep_TwoBitmasks (qs : Qubit[], bits1 : Int[], bits2 : Int[]) : Unit
    is Adj {
        
        let firstDiff = FindFirstDiff_Reference(bits1, bits2);
        H(qs[firstDiff]);
            
        for (i in 0 .. Length(qs) - 1) {
            if (bits1[i] == bits2[i]) {
                if (bits1[i] == 1) {
                    X(qs[i]);
                }
            } else {
                if (i > firstDiff) {
                    CNOT(qs[firstDiff], qs[i]);
                    if (bits1[i] != bits1[firstDiff]) {
                        X(qs[i]);
                    }
                }
            }
        }
    }
    
    operation TestParityOnState (statePrep : (Qubit[] => Unit is Adj), parity : Int, stateStr : String) : Unit {
        
        using (register = Qubit[9]) {
            // prepare basis state to test on
            statePrep(register);

            ResetOracleCallsCount();

            let res = MeasureParity(register);
            
            // check that the returned parity is correct
            Fact((res == Zero) == (parity == 0), $"Failed on {stateStr}.");
            
            // check that the state has not been modified
            Adjoint statePrep(register);
            AssertAllZero(register);

            let nm = GetOracleCallsCount(M) + GetOracleCallsCount(Measure);
            Fact(nm <= 1, $"You are allowed to do at most one measurement, and you did {nm}");
        }
    }
    
    operation P1_Measurement () : Unit {
        // test on all basis states
        for (bits in 0 .. 9) {
            let bitsStr = ToString_Bitmask(bits);
            TestParityOnState(StatePrep_Bitmask(_, bits), Parity(bits), $"basis state |{bitsStr}⟩");
        }
        
        // test on all superpositions of two basis states of the same parity
        for (b1 in 0 .. 9) {
            let bits1 = IntToBoolArray(b1);
            let bitsStr1 = ToString_Bitmask(b1);
            
            for (b2 in b1 + 1 .. 9) {
                if (Parity(b1) == Parity(b2)) {
                    let bits2 = IntToBoolArray(b2);
                    let bitsStr2 = ToString_Bitmask(b2);
                    let p = Parity(b1);
                    Message($"Testing on |{bitsStr1}⟩ + |{bitsStr2}⟩ with parity {p}");
                    TestParityOnState(StatePrep_TwoBitmasks(_, bits1, bits2), Parity(b1), $"state |{bitsStr1}⟩ + |{bitsStr2}⟩");
                }
            }
        }
    }
    
    operation AssertEqualOnZeroState (
        statePrep : (Qubit[] => Unit is Adj), 
        testImpl : (Qubit[] => Unit), 
        refImpl : (Qubit[] => Unit is Adj)) : Unit {
        using (qs = Qubit[3]) {
            // prepare state
            statePrep(qs);
            
            // apply operation that needs to be tested
            testImpl(qs);
            
            // apply adjoint reference operation and adjoint state prep
            Adjoint refImpl(qs);
            Adjoint statePrep(qs);
            
            // assert that all qubits end up in |0⟩ state
            AssertAllZero(qs);
        }
    }
    
    
    operation StatePrep_Rotate (qs : Qubit[], alpha : Double) : Unit
    is Adj {        
        Ry(2.0 * alpha, qs[0]);
    }
    
    operation P2_Encoding () : Unit {
        for (i in 0 .. 810) {
            let alpha = ((2.0 * PI()) * IntAsDouble(i)) / 810.0;
            AssertEqualOnZeroState(StatePrep_Rotate(_, alpha), Encode, Encode_Reference);
        }
    }
    
    operation StatePrep_WithError (qs : Qubit[], alpha : Double, hasError : Bool) : Unit
    is Adj {
        
        StatePrep_Rotate(qs, alpha);
        Encode_Reference(qs);
            
        if (hasError) {
            X(qs[0]);
        }
    }
    
    operation P3_QECC () : Unit {
        using (register = Qubit[9]) {
            for (i in 0 .. 810) {
                let alpha = ((2.0 * PI()) * IntAsDouble(i)) / 810.0;
                StatePrep_WithError(register, alpha, false);
                EqualityFactR(DetectErrorOnLeftQubit(register), Zero, "Error on X");
                Adjoint StatePrep_WithError(register, alpha, false);
                AssertAllZero(register);
                StatePrep_WithError(register, alpha, true);
                EqualityFactR(DetectErrorOnLeftQubit(register), One, "Error on X");
                Adjoint StatePrep_WithError(register, alpha, true);
                AssertAllZero(register);
            }
        }
    }
    
    operation BindErrorCorrectionRoundImpl (
        encoder : (Qubit[] => Unit is Adj), 
        error : Pauli[], 
        logicalOp : (Qubit[] => Unit), 
        correction : (Qubit[] => Unit), 
        dataRegister : Qubit[]) : Unit {
        
        using (auxiliary = Qubit[2]) {
            let register = dataRegister + auxiliary;
            
            // encode the logical qubit (dataRegister) into physical representation (register)
            encoder(register);
            
            // apply error (or no error)
            ApplyPauli(error, register);
            
            // perform logical operation on (possibly erroneous) state
            logicalOp(register);
            
            // apply correction to get the state back to correct one
            correction(register);
            
            // apply decoding to get back to 1-qubit state
            Adjoint encoder(register);
            AssertAllZero(auxiliary);
        }
    }
    
    
    function BindErrorCorrectionRound (
        encoder : (Qubit[] => Unit is Adj), 
        error : Pauli[], 
        logicalOp : (Qubit[] => Unit), 
        correction : (Qubit[] => Unit)) : (Qubit[] => Unit) {
        
        return BindErrorCorrectionRoundImpl(encoder, error, logicalOp, correction, _);
    }
    
    
    // list of errors which can be corrected by the code (the first element corresponds to no error)
    function PauliErrors () : Pauli[][] {
        return [[PauliI, PauliI, PauliI,PauliI,PauliI,PauliI,PauliI,PauliI,PauliI], 
                [PauliX, PauliI, PauliI,PauliI,PauliI,PauliI,PauliI,PauliI,PauliI], 
                [PauliI, PauliX, PauliI,PauliI,PauliI,PauliI,PauliI,PauliI,PauliI], 
                [PauliI, PauliI, PauliX,PauliI,PauliI,PauliI,PauliI,PauliI,PauliI],
                [PauliI, PauliI, PauliI,PauliX,PauliI,PauliI,PauliI,PauliI,PauliI],
                [PauliI, PauliI, PauliI,PauliI,PauliX,PauliI,PauliI,PauliI,PauliI],
                [PauliI, PauliI, PauliI,PauliI,PauliI,PauliX,PauliI,PauliI,PauliI],
                [PauliI, PauliI, PauliI,PauliI,PauliI,PauliI,PauliX,PauliI,PauliI],
                [PauliI, PauliI, PauliI,PauliI,PauliI,PauliI,PauliI,PauliX,PauliI],
                [PauliI, PauliI, PauliI,PauliI,PauliI,PauliI,PauliI,PauliI,PauliX]];
    }
    
    
    operation P4_QECC_LeftQubit () : Unit {
        let partialBind = BindErrorCorrectionRound(Encode_Reference, _, NoOp<Qubit[]>, CorrectErrorOnLeftQubit);
        let errors = PauliErrors();
        
        for (idxError in 0 .. 1) {
            AssertOperationsEqualReferenced(1, partialBind(errors[idxError]), NoOp<Qubit[]>);
        }
    }
    operation P5_QECC_AnyQubit () : Unit {
        let errors = PauliErrors();
        
        using (register = Qubit[9]) {
            
            for (idxError in 0 .. Length(errors) - 1) {
                let θ = RandomReal(12);
                let statePrep = BoundCA([H, Rz(θ, _)]);
                mutable errorStr = "no error";
                if (idxError > 0) {
                    set errorStr = $"error on qubit {idxError}";
                }
                
                Message($"Testing with {errorStr}.");
                statePrep(Head(register));
                Encode_Reference(register);
                ApplyPauli(errors[idxError], register);
                EqualityFactI(DetectErrorOnAnyQubit(register), idxError, $"Failed on state with {errorStr}.");
                ApplyPauli(errors[idxError], register);
                Adjoint Encode_Reference(register);
                Adjoint statePrep(Head(register));
                AssertAllZero(register);
            }
        }
    }
     operation P6_CorrectErrorOnAnyQubit () : Unit {
        
        let partialBind = BindErrorCorrectionRound(Encode_Reference, _, NoOp<Qubit[]>, CorrectErrorOnAnyQubit);
        let errors = PauliErrors();
        
        for (idxError in 0 .. Length(errors) - 1) {
            Message($"{errors[idxError]}");
            AssertOperationsEqualReferenced(1, partialBind(errors[idxError]), NoOp<Qubit[]>);
        }
    }
    operation P7_LogicalX () : Unit {
        
        let partialBind = BindErrorCorrectionRound(Encode_Reference, _, LogicalX, CorrectErrorOnAnyQubit_Reference);
        let errors = PauliErrors();
        
        for (idxError in 0 .. Length(errors) - 1) {
            Message($"{errors[idxError]}");
            AssertOperationsEqualReferenced(1, partialBind(errors[idxError]), ApplyPauli([PauliX], _));
        }
    }
    operation P8_LogicalZ () : Unit {
        
        let partialBind = BindErrorCorrectionRound(Encode_Reference, _, LogicalZ, CorrectErrorOnAnyQubit_Reference);
        let errors = PauliErrors();
        
        for (idxError in 0 .. Length(errors) - 1) {
            Message($"{errors[idxError]}");
            AssertOperationsEqualReferenced(1, partialBind(errors[idxError]), ApplyToEachA(Z, _));
        }
    }
    
}