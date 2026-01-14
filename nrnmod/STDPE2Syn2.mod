: STDP by Hines, changed to dual exponential (BPG 6-1-09)
: Modified by BPG 13-12-08
: Limited weights: max weight is wmax and min weight is wmin
: (initial weight is specified by netconn - usually set to wmin)
: Rhythmic GABAB suppresses conductance and promotes plasticity.
: When GABAB is low, conductance is high and plasticity is off.
:
: Modified by Ieva Kerseviciute 2025-04-10
: Change initial synaptic weight.

NEURON {
	POINT_PROCESS STDPE2bis
	RANGE tau1, tau2, e, i, d, p, dtau, ptau, thresh, wmax, wmin, srcnt1, srcnt2, srcnt3, srcnt4
	RANGE g, gbdel, gblen, gbint, gscale, factor,dshift,dM,dV,B,C
	NONSPECIFIC_CURRENT i

	RANGE initial_weight : The initial synaptic weight
}

UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(uS) = (microsiemens)
}

PARAMETER {
    srcnt1 = 0
    srcnt2 = 0
    srcnt3 = 0
    srcnt4 = 0

    tau1 = 0.5 (ms) <1e-9,1e9>
    tau2 = 3 (ms) <1e-9,1e9>

    e = 0 (mV)

    pi = 3.14159

    wmax = 0.008 (uS)
    wmin = 0 (uS) : not used - use netconn weight instead (BPG)

    d = 0.3 : depression factor (multiplicative to prevent < 0)
    p = 1 : potentiation factor (additive, non-saturating)

    dM = -24 (ms)
    dV = 6.32 (ms)
    ptau = 2 (ms) : Nishiyama2000 10

    thresh = -20 (mV) : postsynaptic voltage threshold

    gbdel = 50 (ms) <1e-9,1e9> : initial GABAB off interval (ms)
    gbint = 125 (ms) <1e-9,1e9> : GABAB off interval (ms)
    gblen = 125 (ms) <1e-9,1e9> : GABAB on length (ms)
    gscale = 1	: relative suppression by GABAB

    dshift = 0 (ms)

    initial_weight = 0 (uS)
}

ASSIGNED {
	v (mV)
	i (nA)
	tpost (ms)
	on
	g (uS)
	gs
	factor
}

STATE {
	C (uS)
	B (uS)
}

INITIAL {
	LOCAL tp
	if (tau1/tau2 > .9999) {
		tau1 = .9999*tau2
	}
	C = 0
	B = 0
	tp = (tau1*tau2)/(tau2 - tau1) * log(tau2/tau1)
	factor = -exp(-tp/tau1) + exp(-tp/tau2)
	factor = 1/factor
	gs=1
	on=1	: initially not plastic
	tpost = -1e9
	net_send(0, 1)
	:net_send(gbdel, 3)	: initial GABAB off period
}

BREAKPOINT {
	SOLVE state METHOD cnexp
	g = B - C
    i = g * gs * (v - e)
}

DERIVATIVE state {
	C' = -C/tau1
	B' = -B/tau2
}

NET_RECEIVE(w (uS), A, tpre (ms)) {
    INITIAL {
        A = initial_weight
        tpre = -1e9
    }

    if (flag == 0) { : Presynaptic spike (after last post so depress)
        C = C + (w + A) * factor
        B = B + (w + A) * factor
        tpre = t

        if (on == 1) {
            A = A * (1 - (d * exp(-((tpost - t) - dM)^2 / (2 * dV * dV))))
        }
    }
    else if (flag == 2 && on == 1) { : Postsynaptic spike
        tpost = t

        FOR_NETCONS(w1, A1, tp) { : Also can hide NET_RECEIVE args
            A1 = A1 + (2 * w1 - w1 - A1) * p * exp((tp - t) / ptau)
        }
    }
    else if (flag == 1) { : flag == 1 from INITIAL block
        WATCH (v > thresh) 2
    }
    else if (flag == 3) { : Plasticity control
        if (on == 0) { : Start plasticity
            on = 1
            gs = gscale
            net_send(gblen, 3)
        }
        else { : End burst
            on = 0
            gs = 1
            net_send(gbint, 3)
        }
    }
}
