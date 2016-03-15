from etalia.consumers.models import ConsumerPubmed, ConsumerArxiv
from etalia.library.models import Journal

journals = ['Nature',
            'Science',
            'The New England journal of medicine']

cp = ConsumerPubmed.objects.first()

# deactivate all journals
cp.deactivate_all()
# activate only those with name in journals
for journal_name in journals:
    js = Journal.objects.filter(title=journal_name)
    for j in js:
        cp.activate_journal(j)
        cp.populate_journal(j.pk)


