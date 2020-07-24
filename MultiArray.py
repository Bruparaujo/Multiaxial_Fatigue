class MultiArray(object):
    def __init__(self, *mains):

        self.count = 0
        self.index = 0

        for main in mains:
            setattr(self, 'array_{}'.format(self.count), main)
            self.count += 1

    # Para obter algum array de dentro do multiArrayInstance
    def get(self, n):
        return getattr(self, 'array_{}'.format(n))

    def first(self):
        return self.get(0)

    def last(self):
        return self.get(self.count - 1)

    # __iter__ + __next__ sao para usar (for var in multiArrayInstance)
    def __iter__(self):
        return self

    def __next__(self):
        if self.index >= self.count:
            self.index = 0
            raise StopIteration
            
        self.index += 1
        return self.get(self.index - 1)
    
    # Descreve como funciona o print(multiArrayInstance)
    def __str__(self):
        ret = ''
        for arr in self:
            ret = '{}{}'.format(ret, '{} / {} - {}\n'.format(self.index, self.count, arr))
        return ret

# Self test
if __name__ == '__main__':
    # 1
    # Arrange
    entryTest11 = [1, 2, 3, 4, 5]
    # Act
    testA = MultiArray(entryTest11)
    # Assert
    print(testA.get(0))
    print(testA.count)

    # 2
    entryTest12 = [5, 4, 3, 2, 1]
    testB = MultiArray(entryTest11, entryTest12)
    print(testB.get(0))
    print(testB.get(1))
    print(testB.get(0)+testB.get(1))
    asd
    # 3
    for array in testB:
        print('{} / {} - {}'.format(testB.index, testB.count, array))

    # 4
    auxs = MultiArray( *([] for _ in range(4) ) )
    for i, aux in enumerate(auxs):
        aux.append(i)
        print(aux)
    print(auxs)
    
    # 5
    auxs.get(2).append('a')
    print(auxs.get(2))

    # 6
    print(auxs.first())

    # 7
    print(auxs.last())

