define([
    'jquery',
    'app/ui/api',
    'app/util/utils',
    'app/ui/layout',
    'app/ui/controls',
    'app/ui/list',
    'app/ui/detail',
    'endless',
    'bootstrap'
], function($, api, utils, layout, controls, List, Detail) {

    var $body, $list,
        list = new List(),
        detail = new Detail(),
        listType, $clearTrashButton,
        countsHandler, pinHandler, addHandler, trashHandler, restoreHandler;

    pinHandler = function(e, result) {
        $('.thumb[data-paper-id=' + result.getId() + ']')
            .find('.thumb-pin')
            .toggleClass('active', result.isPinned());
    };

    addHandler = restoreHandler = function(e, result) {
        if (result.isAdded()) {
            $('.thumb[data-paper-id=' + result.getId() + ']').remove();
        } else {
            var $button = $(
                '.thumb[data-paper-id=' + result.getId() + '] .thumb-library-add, ' +
                '.thumb[data-paper-id=' + result.getId() + '] .thumb-library-restore '
            ).eq(0);
            toggleLibraryAddOrTrash($button, false);
        }
    };

    trashHandler = function(e, result) {
        if (result.isTrashed()) {
            $('.thumb[data-paper-id=' + result.getId() + ']').remove();
        } else {
            var $button = $('.thumb[data-paper-id=' + result.getId() + '] .thumb-library-trash');
            toggleLibraryAddOrTrash($button, true);
        }
    };

    countsHandler = function(e, result) {
        $('.user-library span').text(result.getLibraryCount());
        $('.user-library-pins span').text(result.getPinCount());
        $('.user-library-trash span').text(result.getTrashCount());

        if (0 < result.getTrashCount()) {
            $clearTrashButton.show();
        } else {
            $clearTrashButton.hide();
        }
    };

    function toggleLibraryAddOrTrash($button, added) {
        if (added) {
            $button
                .removeClass('thumb-library-add thumb-library-restore')
                .addClass('thumb-library-trash');

            utils.restoreLoadingButton($button, 'eai-library-trash');
        } else {
            $button
                .removeClass('thumb-library-trash')
                .addClass('thumb-library-restore');

            utils.restoreLoadingButton($button, 'eai-library-add');
        }
    }

    $(function() {

        $body = $('body');
        $list = $('#list');

        list.init();
        detail.init();

        listType = $list.data('type');
        $clearTrashButton = $('#user-library-trash-clear');

        if (listType == 'pin') {
            pinHandler = function(e, result) {
                if (!result.isPinned()) {
                    $('.thumb[data-paper-id=' + result.getId() + ']').remove();
                }
            };
            addHandler = restoreHandler = function(e, result) {
                var $button = $(
                    '.thumb[data-paper-id=' + result.getId() + '] .thumb-library-add, ' +
                    '.thumb[data-paper-id=' + result.getId() + '] .thumb-library-restore '
                ).eq(0);
                if (result.isAdded()) {
                    toggleLibraryAddOrTrash($button, true);
                } else {
                    toggleLibraryAddOrTrash($button, false);
                }
            };
            trashHandler = function(e, result) {
                var $button = $('.thumb[data-paper-id=' + result.getId() + '] .thumb-library-trash');
                if (result.isTrashed()) {
                    toggleLibraryAddOrTrash($button, false);
                } else {
                    toggleLibraryAddOrTrash($button, true);
                }
            };
        }

        $clearTrashButton.on('click', function(e) {
            if (confirm("Are you sure you want to clear the trash ?")) {
                api.clearTrash();
            }

            e.preventDefault();
            return false;
        });

        $body
            .on(
                'etalia.control.search.change ' +
                'etalia.control.timespan.change ' +
                'etalia.control.cluster.change ' +
                'etalia.control.pinned.change ' +
                'etalia.control.filters.change',
                function() {
                    list.load(controls.getStates());
                })
            .on('etalia.list.load', function(e, data) {
                if (data.hasOwnProperty('controlsStates')) {
                    controls.setStates(data.controlsStates);
                }
            })

            .on('etalia.publication.pin', pinHandler)
            .on('etalia.publication.add', addHandler)
            .on('etalia.publication.trash', trashHandler)
            .on('etalia.publication.restore', restoreHandler)

            .on('etalia.publication.pin ' +
                'etalia.publication.add ' +
                'etalia.publication.trash ' +
                'etalia.publication.restore ' +
                'etalia.publication.trash-clear',
                countsHandler)

            .on('etalia.publication.trash-clear', function() {
                $('.thumb-list').empty();
                $clearTrashButton.hide();
            })

            .on('etalia.detail.loading', function() {
                layout.setBusy();
            })
            .on('etalia.detail.loaded', function() {
                layout.setAvailable();
            });

        $list.on('click', '.thumb .title a', function(e) {
            detail.load($(e.target).closest('.thumb'));

            e.preventDefault();
            return false;
        });

        // Endless scroll
        $.endlessPaginate({
            paginateOnScroll: true,
            paginateOnScrollMargin: 10,
            onClick: function (context) {
                context.extraData = JSON.stringify(controls.getStates());
                context.urlOrigin = window.location.href;
            }
        });
    });
});
