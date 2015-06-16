import logging
from config.celery import celery_app as app

logger = logging.getLogger(__name__)

# This is very non ideal, but I cannot make celery call the backend class method
# because serialization failed during in kombu module.
#
@app.task(name='mendeley.update_lib')
def update_lib_mendeley_static(user, mendeley_session, name, parser,
                               get_or_create_entry, associate_paper,
                               associate_journal, CHUNK_SIZE,
                               create_lib_stats):

    # Init
    new = True
    count = 0
    user.lib.set_lib_syncing()

    # retrieve list of documents per page
    page = mendeley_session.documents.list(
        page_size=CHUNK_SIZE,
        sort='created',
        order='desc',
        view='all')

    # loop through documents/page until end or find already registered doc
    while True:
        for item in page.items:
            entry = parser.parse(item)
            paper, journal = get_or_create_entry(entry)

            if paper:
                logger.info(
                    '+ Entry: {ids} from {user} / {backend}'.format(
                        ids=paper.print_ids,
                        user=user.email,
                        backend=name))
                new = associate_paper(paper, user, entry['user_info']) and new
                if new:
                    count += 1
                else:
                    # escape when reaching already uploaded references
                    break
                if journal:
                    associate_journal(journal, user)
            else:
                logger.info(
                    '- Item: {type_} from {user} / {backend}'.format(
                        type_=item.type,
                        user=user.email,
                        backend=name))
        if page.next_page:
            page = page.next_page
        else:
            break

    # update UserLib and Stats
    create_lib_stats(user, count)
    user.lib.set_lib_idle()
    return count
