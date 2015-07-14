from nlp.models import Model
from nlp.utils import TaggedDocumentsIterator
from time import time

t0 = time()

model = Model.objects.create(name='dbow',
                             dm=0,
                             size=100,
                             negative=5,
                             hs=0,
                             min_count=2,
                             workers=3)

docs = TaggedDocumentsIterator('nlp/data/')
documents = [doc for doc in docs]
model.build_vocab(documents)
model.train(documents, passes=10)
t1 = time()
print('Training in {0}'.format(t1 - t0))
model.save()
model.populate_library()
print('Saving in {0}'.format(time() - t1))
