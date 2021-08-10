

def test_exp_add_to_avoid_zero():
    '''
    '''
    exp_df = pd.DataFrame(index=[0,1,2,3],columns=['a','b'])

    # ...

def test_find_min_after_zero():
    '''
    '''
    df = pd.DataFrame({'a':[1,2,3,4],'b':[0,0,1,3]})
    res = find_min_after_zero(df)
    assert(res==1)