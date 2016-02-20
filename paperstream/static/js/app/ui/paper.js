define(['jquery', 'app/api', 'app/util/utils', 'app/ui/controls', 'app/ui/list', 'app/util/sticky'], function ($, api, utils, controls, List, Sticky) {

    var Paper = function (options) {

        this.config = $.extend({
            debug: false,
            element: '#detail',
            document: '.document',
            actions: '.inner > .actions',
            neighbors: '.neighbors-thumbs'
        }, options);

        this.$element = $(this.config.element);
        this.$document = this.$element.find(this.config.document);

        this.id = null;
        this.$actions = null;
        this.actions = null;

        this.neighbors = null;
    };

    function toggleLibraryAddOrTrash($button, added) {
        if (added) {
            $button
                .removeClass('detail-library-add detail-library-restore')
                .addClass('detail-library-trash');

            utils.restoreLoadingButton($button, 'eai-library-trash');
        } else {
            $button
                .removeClass('detail-library-trash')
                .addClass('detail-library-restore');

            utils.restoreLoadingButton($button, 'eai-library-add');
        }
    }

    Paper.prototype.log = function() {
        if (this.config.debug) {
            console.log('[Paper] ' + arguments[0], Array.prototype.splice.call(arguments, 1));
        }
    };

    Paper.prototype.getId = function() {
        return this.id;
    };

    Paper.prototype.init = function() {
        this.log('init');

        // API events
        var that = this;
        $('body')
            .on('etalia.publication.pin', function(e, result) {
                if (that.id == result.getId()) {
                    that.$actions
                        .find('.detail-pin')
                        .toggleClass('active', result.isPinned());
                }
            })
            .on('etalia.publication.add', function(e, result) {
                if (that.id == result.getId() && result.isAdded()) {
                    var $button = that.$actions.find('.detail-library-add');
                    toggleLibraryAddOrTrash($button, true);
                }
            })
            .on('etalia.publication.trash', function(e, result) {
                if (that.id == result.getId() && result.isTrashed()) {
                    var $button = that.$actions.find('.detail-library-trash');
                    toggleLibraryAddOrTrash($button, false);
                }
            })
            .on('etalia.publication.restore', function(e, result) {
                if (that.id == result.getId() && result.isAdded()) {
                    var $button = that.$actions.find('.detail-library-restore');
                    toggleLibraryAddOrTrash($button, true);
                }
            });

        this.$document
            .on('click', '.neighbors-timespan-selector a', function(e) {
                var $timespan = $(e.target).closest('a');

                $timespan.closest('.neighbors-timespan-selector').find('li').removeClass('active');
                $timespan.closest('li').addClass('active');

                if (that.neighbors) {
                    that.neighbors.load({'time_span': $timespan.data('value')});
                }

                e.preventDefault();
                return false;
            })
            .on('click', '.detail-questions', function(e) {
                var top = that.$document.find('.questions-thumbs').position().top - 15;
                that.$document.scrollTop(top);

                e.preventDefault();
                return false;
            })
            .on('click', '.detail-neighbors', function(e) {
                var top = that.$document.find('.neighbors-thumbs').position().top - 15;
                that.$document.scrollTop(top);

                e.preventDefault();
                return false;
            })
            .on('click', '.detail-pin', function(e) {
                if (!that.id) throw 'Undefined paper id';

                api.pin(that.id);

                e.preventDefault();
                return false;
            })
            .on('click', '.detail-ban', function(e) {
                if (!that.id) throw 'Undefined paper id';

                api.ban(that.id);

                e.preventDefault();
                return false;
            })
            .on('click', '.detail-library-add:visible', function(e) {
                if (!that.id) throw 'Undefined paper id';

                var $button = $(e.target).closest('.detail-library-add');
                api.add(that.id, function() {
                    utils.restoreLoadingButton($button, 'eai-library-add');
                });

                e.preventDefault();
                return false;
            })
            .on('click', '.detail-library-trash:visible', function(e) {
                if (!that.id) throw 'Undefined paper id';

                var $button = $(e.target).closest('.detail-library-trash');
                api.trash(that.id, function() {
                    utils.restoreLoadingButton($button, 'eai-library-trash');
                });

                e.preventDefault();
                return false;
            })
            .on('click', '.detail-library-restore:visible', function(e) {
                if (!that.id) throw 'Undefined paper id';

                var $button = $(e.target).closest('.detail-library-restore');
                api.restore(that.id, function() {
                    utils.restoreLoadingButton($button, 'eai-library-add');
                });

                e.preventDefault();
                return false;
            })
            .on('click', '.detail-twitter:visible', function(e) {
                var $actions = $('#detail-actions'),
                    url = 'https://twitter.com/intent/tweet/'
                        + '?text=' + $actions.data('paper-title')
                        + '&url=' + $actions.data('paper-url')
                        + '&via=etalia';
                        //+ '&hashtags=web,development';

                utils.popup(url, 'share-popup');

                e.preventDefault();
                return false;
            })
            .on('click', '.detail-google-plus:visible', function(e) {
                var $actions = $('#detail-actions'),
                    url = 'https://plus.google.com/share'
                        + '?url=' + $actions.data('paper-url');

                utils.popup(url, 'share-popup');

                e.preventDefault();
                return false;
            });

        return this;
    };

    Paper.prototype.update = function() {
        this.log('update');

        this.$actions = this.$element.find(this.config.actions);
        utils.bindLoadingButtons(this.$actions);

        this.id = parseInt(this.$actions.data('paper-id'));

        // Sticky action
        this.actions = new Sticky({
            debug: this.config.debug,
            element: this.$actions,
            parent: this.$document,
            top: 20,
            bottom: 20
        });
        this.actions.enable();

        var timespan = $('.neighbors-timespan-selector li.active a').data('value') || 30;

        this.neighbors = new List({
                debug: this.config.debug,
                element: this.config.element + ' ' + this.config.neighbors
            })
            .init()
            .load({'time_span': timespan});

        return this;
    };

    Paper.prototype.clear = function() {
        if (this.actions) {
            this.actions.disable();
            this.actions = null;
        }

        if (this.$actions) {
            this.$actions = null;
        }

        this.neighbors = null;

        if (this.id) {
            this.id = null;
        }
    };

    return Paper;
});
