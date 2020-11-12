import numpy as np

def over(xx, yy):
    return np.sum(xx * yy)

#def err(xx, yy):
#    return abs(over(xx, yy)) ** 0.5/ np.sum(yy**2)**0.5

def err(xx, yy):
    return np.sum((xx - yy)**2) ** 0.5/ np.sum(yy**2)**0.5


def shape_err(xx, yy):
    xx1 = xx / over(xx, xx)**0.5
    yy1 = yy / over(yy, yy)**0.5

    return abs(over(xx1, yy1))


def err_mult(xx, yy):
    return abs(np.sum((xx * yy))) ** 0.5/ np.sum(yy**2)**0.5
    #xx1 = xx / over(xx, xx)**0.5
    #yy1 = yy / over(yy, yy)**0.5
    #return abs(over(xx1, yy1))   
def over(xx, yy):
    return np.sum(xx * yy)

def create_overlap_mat(*states4 ):
    #if states2.all() ==  None:
    #import ipdb
    #ipdb.set_trace()
    if len(states4) == 1:
        states = states4[0]
        n_states = states.shape[0]
        over_mat = np.zeros((n_states,n_states))
        for i_state in range(n_states):
            for j_state in range(n_states):
                over_mat[i_state,j_state] =  over(states[i_state], states[j_state])
        return over_mat
    else:
        #import ipdb
        #ipdb.set_trace()
        states = states4[0]
        states2 = states4[1]
        #states2  = np.array(states3)
        n_states = states.shape[0]
        m_states = states2.shape[0]
        over_mat = np.zeros((n_states,m_states))
        for i_state in range(n_states):
            for j_state in range(m_states):
                over_mat[i_state,j_state] =  over(states[i_state], states2[j_state])
        return over_mat

