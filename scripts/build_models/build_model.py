import yaml
from nlp.models import Model, LSH
from library.models import Paper
from core.constants import TIME_LAPSE_CHOICES

path_to_model_file = 'scripts/build_models/models.yaml'
with open(path_to_model_file) as stream:
    MODELS = yaml.load(stream)

def build(model_name):

    check_unique_name(MODELS)
    model_args = [model for model in MODELS if model['name'] == model_name][0]

    # Initiate model
    model = Model.objects.create(**model_args)
    # dump papers data
    papers = Paper.objects.all()
    model.dump(papers)
    # Build vocabulary and phraser
    docs = model.load_documents()
    phraser = model.build_phraser(docs, min_count=2)
    docs = model.load_documents(phraser=phraser)
    model.build_vocab(docs)
    # Train, save and set_active
    model.train(docs, passes=10, shuffle_=True)
    model.save()
    model.set_active()
    # Populate library
    model.save_journal_vec_from_bulk()
    model.save_paper_vec_from_bulk()
    # Build related LSH
    for (time_lapse, _) in TIME_LAPSE_CHOICES:
        LSH.objects.create(model=model,
                           time_lapse=time_lapse)


def check_unique_name(models):
    names = [(i, models['name']) for i, model in enumerate(models)]
    if not len(names) == len(set(names)):
        raise AssertionError('MODELS has duplicated names')