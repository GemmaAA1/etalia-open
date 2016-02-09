define(['jquery', 'app/ui/list', 'bootstrap'], function($) {

    var listType, countsHandler, pinHandler, addHandler, trashHandler, restoreHandler;

    pinHandler = function(e, data) {
        if (data.hasOwnProperty('is_liked')) {
            $('.thumb[data-id=' + data.id + ']')
                .find('.thumb-pin')
                .toggleClass('active', data['is_liked']);
        }
    };

    addHandler = trashHandler = restoreHandler = function(e, data) {
        if (data.hasOwnProperty('success') && data['success']) {
            $('.thumb[data-id=' + data.id + ']').remove();
        }
    };

    countsHandler = function(e, data) {
        if (data.hasOwnProperty('library_counter')) {
            $('.user-library span').text(data['library_counter']);
        }
        if (data.hasOwnProperty('likes_counter')) {
            $('.user-library-pins span').text(data['likes_counter']);
        }
        if (data.hasOwnProperty('trash_counter')) {
            $('.user-library-trash span').text(data['trash_counter']);
        }
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
            pinHandler = function(e, data) {
                if (data.hasOwnProperty('is_liked') && !data['is_liked']) {
                    $('.thumb[data-id=' + data.id + ']').remove();
                }
            };
            addHandler = function(e, data) {
                if (data.hasOwnProperty('success') && data['success']) {
                    var $button = $('.thumb[data-id=' + data.id + '] .thumb-library-add');
                    toggleLibraryAddOrTrash($button, true);
                }
            };
            trashHandler = function(e, data) {
                if (data.hasOwnProperty('success') && data['success']) {
                    var $button = $('.thumb[data-id=' + data.id + '] .thumb-library-trash');
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
