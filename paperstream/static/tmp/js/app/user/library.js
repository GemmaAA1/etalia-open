define(['jquery', 'app/ui/list', 'bootstrap'], function($) {

    var listType, countsHandler, pinHandler, addHandler, trashHandler, restoreHandler;

    pinHandler = function(e, result) {
        $('.thumb[data-id=' + result.getId() + ']')
            .find('.thumb-pin')
            .toggleClass('active', result.isPinned());
    };

    addHandler = restoreHandler = function(e, result) {
        if (result.isAdded()) {
            $('.thumb[data-id=' + result.getId() + ']').remove();
        }
    };

    trashHandler = function(e, result) {
        if (result.isTrashed()) {
            $('.thumb[data-id=' + result.getId() + ']').remove();
        }
    };

    countsHandler = function(e, result) {
        $('.user-library span').text(result.getLibraryCount());
        $('.user-library-pins span').text(result.getPinCount());
        $('.user-library-trash span').text(result.getTrashCount());
    };

    function toggleLibraryAddOrTrash($button, added) {
        if (added) {
            $button
                .removeClass('thumb-library-add')
                .addClass('thumb-library-trash')
                .find('.eai')
                    .removeClass('eai-library-add')
                    .addClass('eai-library-trash');
        } else {
            $button
                .removeClass('thumb-library-trash')
                .addClass('thumb-library-restore')
                .find('.eai')
                    .removeClass('eai-library-trash')
                    .addClass('eai-library-add');
        }
    }

    $(function() {

        listType = $('#list').data('type');

        if (listType == 'pin') {
            pinHandler = function(e, result) {
                if (!result.isPinned()) {
                    $('.thumb[data-id=' + result.getId() + ']').remove();
                }
            };
            addHandler = function(e, result) {
                if (result.isAdded()) {
                    var $button = $('.thumb[data-id=' + result.getId() + '] .thumb-library-add');
                    toggleLibraryAddOrTrash($button, true);
                }
            };
            trashHandler = function(e, result) {
                if (result.isTrashed()) {
                    var $button = $('.thumb[data-id=' + result.getId() + '] .thumb-library-trash');
                    toggleLibraryAddOrTrash($button, false);
                }
            };
        }

        $('body')
            .on('etalia.publication.pin', pinHandler)
            .on('etalia.publication.add', addHandler)
            .on('etalia.publication.trash', trashHandler)
            .on('etalia.publication.restore', restoreHandler)

            .on('etalia.publication.pin ' +
                'etalia.publication.add ' +
                'etalia.publication.trash ' +
                'etalia.publication.restore',
                countsHandler);
    });
});
