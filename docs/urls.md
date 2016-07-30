## URL endpoints

    /api/v1/	rest_framework.routers.APIRoot	api:api-root	
    /api/v1/.<format>/	rest_framework.routers.APIRoot	api:api-root	
    /api/v1/library/authors.<format>/	etalia.library.api.views.AuthorViewSet	api:author-list	
    /api/v1/library/authors/	etalia.library.api.views.AuthorViewSet	api:author-list	
    /api/v1/library/authors/<pk>.<format>/	etalia.library.api.views.AuthorViewSet	api:author-detail	
    /api/v1/library/authors/<pk>/	etalia.library.api.views.AuthorViewSet	api:author-detail	
    /api/v1/library/journals.<format>/	etalia.library.api.views.JournalViewSet	api:journal-list	
    /api/v1/library/journals/	etalia.library.api.views.JournalViewSet	api:journal-list	
    /api/v1/library/journals/<pk>.<format>/	etalia.library.api.views.JournalViewSet	api:journal-detail	
    /api/v1/library/journals/<pk>/	etalia.library.api.views.JournalViewSet	api:journal-detail	
    /api/v1/library/papers.<format>/	etalia.library.api.views.PaperViewSet	api:paper-list	
    /api/v1/library/papers/	etalia.library.api.views.PaperViewSet	api:paper-list	
    /api/v1/library/papers/<pk>.<format>/	etalia.library.api.views.PaperViewSet	api:paper-detail	
    /api/v1/library/papers/<pk>/	etalia.library.api.views.PaperViewSet	api:paper-detail	
    /api/v1/library/papers/<pk>/neighbors.<format>/	etalia.library.api.views.PaperViewSet	api:paper-neighbors	
    /api/v1/library/papers/<pk>/neighbors/	etalia.library.api.views.PaperViewSet	api:paper-neighbors	
    /api/v1/library/papers/<pk>/related-threads.<format>/	etalia.library.api.views.PaperViewSet	api:paper-related-threads	
    /api/v1/library/papers/<pk>/related-threads/	etalia.library.api.views.PaperViewSet	api:paper-related-threads	
    /api/v1/library/papers/filters.<format>/	etalia.library.api.views.PaperViewSet	api:paper-filters	
    /api/v1/library/papers/filters/	etalia.library.api.views.PaperViewSet	api:paper-filters	
    /api/v1/library/states.<format>/	etalia.library.api.views.PaperStateViewSet	api:paperuser-list	
    /api/v1/library/states/	etalia.library.api.views.PaperStateViewSet	api:paperuser-list	
    /api/v1/library/states/<pk>.<format>/	etalia.library.api.views.PaperStateViewSet	api:paperuser-detail	
    /api/v1/library/states/<pk>/	etalia.library.api.views.PaperStateViewSet	api:paperuser-detail	
    /api/v1/library/states/empty-trash.<format>/	etalia.library.api.views.PaperStateViewSet	api:paperuser-empty-trash	
    /api/v1/library/states/empty-trash/	etalia.library.api.views.PaperStateViewSet	api:paperuser-empty-trash	
    /api/v1/popover/popovers.<format>/	etalia.popovers.api.views.PopOverViewSet	api:popover-list	
    /api/v1/popover/popovers/	etalia.popovers.api.views.PopOverViewSet	api:popover-list	
    /api/v1/popover/popovers/<pk>.<format>/	etalia.popovers.api.views.PopOverViewSet	api:popover-detail	
    /api/v1/popover/popovers/<pk>/	etalia.popovers.api.views.PopOverViewSet	api:popover-detail	
    /api/v1/popover/states.<format>/	etalia.popovers.api.views.PopOverStateViewSet	api:userpopover-list	
    /api/v1/popover/states/	etalia.popovers.api.views.PopOverStateViewSet	api:userpopover-list	
    /api/v1/popover/states/<pk>.<format>/	etalia.popovers.api.views.PopOverStateViewSet	api:userpopover-detail	
    /api/v1/popover/states/<pk>/	etalia.popovers.api.views.PopOverStateViewSet	api:userpopover-detail	
    /api/v1/thread/comments.<format>/	etalia.threads.api.views.ThreadCommentViewSet	api:threadcomment-list	
    /api/v1/thread/comments/	etalia.threads.api.views.ThreadCommentViewSet	api:threadcomment-list	
    /api/v1/thread/comments/<pk>.<format>/	etalia.threads.api.views.ThreadCommentViewSet	api:threadcomment-detail	
    /api/v1/thread/comments/<pk>/	etalia.threads.api.views.ThreadCommentViewSet	api:threadcomment-detail	
    /api/v1/thread/invites.<format>/	etalia.threads.api.views.ThreadUserInviteViewSet	api:threaduserinvite-list	
    /api/v1/thread/invites/	etalia.threads.api.views.ThreadUserInviteViewSet	api:threaduserinvite-list	
    /api/v1/thread/invites/<pk>.<format>/	etalia.threads.api.views.ThreadUserInviteViewSet	api:threaduserinvite-detail	
    /api/v1/thread/invites/<pk>/	etalia.threads.api.views.ThreadUserInviteViewSet	api:threaduserinvite-detail	
    /api/v1/thread/posts.<format>/	etalia.threads.api.views.ThreadPostViewSet	api:threadpost-list	
    /api/v1/thread/posts/	etalia.threads.api.views.ThreadPostViewSet	api:threadpost-list	
    /api/v1/thread/posts/<pk>.<format>/	etalia.threads.api.views.ThreadPostViewSet	api:threadpost-detail	
    /api/v1/thread/posts/<pk>/	etalia.threads.api.views.ThreadPostViewSet	api:threadpost-detail	
    /api/v1/thread/states.<format>/	etalia.threads.api.views.ThreadUserViewSet	api:threaduser-list	
    /api/v1/thread/states/	etalia.threads.api.views.ThreadUserViewSet	api:threaduser-list	
    /api/v1/thread/states/<pk>.<format>/	etalia.threads.api.views.ThreadUserViewSet	api:threaduser-detail	
    /api/v1/thread/states/<pk>/	etalia.threads.api.views.ThreadUserViewSet	api:threaduser-detail	
    /api/v1/thread/threads.<format>/	etalia.threads.api.views.ThreadViewSet	api:thread-list	
    /api/v1/thread/threads/	etalia.threads.api.views.ThreadViewSet	api:thread-list	
    /api/v1/thread/threads/<pk>.<format>/	etalia.threads.api.views.ThreadViewSet	api:thread-detail	
    /api/v1/thread/threads/<pk>/	etalia.threads.api.views.ThreadViewSet	api:thread-detail	
    /api/v1/thread/threads/<pk>/neighbors.<format>/	etalia.threads.api.views.ThreadViewSet	api:thread-neighbors	
    /api/v1/thread/threads/<pk>/neighbors/	etalia.threads.api.views.ThreadViewSet	api:thread-neighbors	
    /api/v1/thread/threads/<pk>/publish.<format>/	etalia.threads.api.views.ThreadViewSet	api:thread-publish	
    /api/v1/thread/threads/<pk>/publish/	etalia.threads.api.views.ThreadViewSet	api:thread-publish	
    /api/v1/thread/threads/filters.<format>/	etalia.threads.api.views.ThreadViewSet	api:thread-filters	
    /api/v1/thread/threads/filters/	etalia.threads.api.views.ThreadViewSet	api:thread-filters	
    /api/v1/user/affiliations.<format>/	etalia.users.api.views.AffiliationViewSet	api:affiliation-list	
    /api/v1/user/affiliations/	etalia.users.api.views.AffiliationViewSet	api:affiliation-list	
    /api/v1/user/affiliations/<pk>.<format>/	etalia.users.api.views.AffiliationViewSet	api:affiliation-detail	
    /api/v1/user/affiliations/<pk>/	etalia.users.api.views.AffiliationViewSet	api:affiliation-detail	
    /api/v1/user/relationships.<format>/	etalia.users.api.views.RelationshipViewSet	api:relationship-list	
    /api/v1/user/relationships/	etalia.users.api.views.RelationshipViewSet	api:relationship-list	
    /api/v1/user/relationships/<pk>.<format>/	etalia.users.api.views.RelationshipViewSet	api:relationship-detail	
    /api/v1/user/relationships/<pk>/	etalia.users.api.views.RelationshipViewSet	api:relationship-detail	
    /api/v1/user/user-lib-papers.<format>/	etalia.users.api.views.UserLibPaperViewSet	api:userlibpaper-list	
    /api/v1/user/user-lib-papers/	etalia.users.api.views.UserLibPaperViewSet	api:userlibpaper-list	
    /api/v1/user/user-lib-papers/<pk>.<format>/	etalia.users.api.views.UserLibPaperViewSet	api:userlibpaper-detail	
    /api/v1/user/user-lib-papers/<pk>/	etalia.users.api.views.UserLibPaperViewSet	api:userlibpaper-detail	
    /api/v1/user/user-libs.<format>/	etalia.users.api.views.UserLibViewSet	api:userlib-list	
    /api/v1/user/user-libs/	etalia.users.api.views.UserLibViewSet	api:userlib-list	
    /api/v1/user/user-libs/<pk>.<format>/	etalia.users.api.views.UserLibViewSet	api:userlib-detail	
    /api/v1/user/user-libs/<pk>/	etalia.users.api.views.UserLibViewSet	api:userlib-detail	
    /api/v1/user/user-libs/<pk>/papers.<format>/	etalia.users.api.views.UserLibViewSet	api:userlib-papers	
    /api/v1/user/user-libs/<pk>/papers/	etalia.users.api.views.UserLibViewSet	api:userlib-papers	
    /api/v1/user/users.<format>/	etalia.users.api.views.UserViewSet	api:user-list	
    /api/v1/user/users/	etalia.users.api.views.UserViewSet	api:user-list	
    /api/v1/user/users/<pk>.<format>/	etalia.users.api.views.UserViewSet	api:user-detail	
    /api/v1/user/users/<pk>/	etalia.users.api.views.UserViewSet	api:user-detail	
    /api/v1/user/users/<pk>/blocked.<format>/	etalia.users.api.views.UserViewSet	api:user-blocked	
    /api/v1/user/users/<pk>/blocked/	etalia.users.api.views.UserViewSet	api:user-blocked	
    /api/v1/user/users/<pk>/followers.<format>/	etalia.users.api.views.UserViewSet	api:user-followers	
    /api/v1/user/users/<pk>/followers/	etalia.users.api.views.UserViewSet	api:user-followers	
    /api/v1/user/users/<pk>/following.<format>/	etalia.users.api.views.UserViewSet	api:user-following	
    /api/v1/user/users/<pk>/following/	etalia.users.api.views.UserViewSet	api:user-following	
    /feeds/	etalia.feeds.views.my_feeds	feeds:my_feeds	login_required
    /feeds/stream/<name>/reset	etalia.feeds.views.reset_stream_view	feeds:reset-stream	login_required
    /feeds/stream/<name>/update	etalia.feeds.views.update_stream_view	feeds:update-stream	login_required
    /feeds/trend/<name>/reset	etalia.feeds.views.reset_trend_view	feeds:reset-trend	login_required
    /feeds/trend/<name>/update	etalia.feeds.views.update_trend_view	feeds:update-trend	login_required
    /help/	etalia.core.views.help	core:help	
    /papers/	etalia.library.views.my_papers	papers:my_papers	
    /support/	etalia.core.views.support	core:support	
    /terms-of-privacy/	etalia.core.views.terms_privacy	core:terms_privacy	
    /terms-of-use/	etalia.core.views.terms_use	core:terms_use	
    /test-failing-task	etalia.core.views.test_failing_task		
    /threads/	etalia.threads.views.my_threads	threads:my_threads	
    /threads/<pk>/	etalia.threads.views.thread	threads:thread	
    /user/avatar/add/	avatar.views.add	user:avatar_add	login_required
    /user/avatar/change/	avatar.views.change	user:avatar_change	login_required
    /user/avatar/delete/	avatar.views.delete	user:avatar_delete	login_required
    /user/avatar/list/<username>/	avatar.views.avatar_gallery	user:avatar_gallery	login_required
    /user/avatar/list/<username>/<id>/	avatar.views.avatar	user:avatar	login_required
    /user/avatar/render_primary/<user>/<size>/	avatar.views.render_primary	user:avatar_render_primary	login_required
    /user/complete/<backend>/	social.apps.django_app.views.complete	social:complete	
    /user/delete	etalia.users.views.UserDeleteView	user:delete-user	
    /user/disconnect/<backend>/	social.apps.django_app.views.disconnect	social:disconnect	
    /user/disconnect/<backend>/<association_id>/	social.apps.django_app.views.disconnect	social:disconnect_individual	
    /user/done/	etalia.users.views.done	user:done	login_required
    /user/email-sent/	etalia.users.views.validation_sent	user:validation-sent	login_required
    /user/info-affiliation/	etalia.users.views.UserAffiliationSignupView	user:require-affiliation	
    /user/info/	etalia.users.views.UserBasicInfoSignupView	user:require-basic-info	
    /user/login/<backend>/	social.apps.django_app.views.auth	social:begin	
    /user/logout/	etalia.users.views.logout	user:logout	login_required
    /user/profile/	etalia.users.views.ProfileView	user:profile	
    /user/send-invite	etalia.users.views.send_invite	user:send-invite	login_required
    /user/settings/	etalia.users.views.SettingsView	user:settings	
    /user/signin/	etalia.users.views.UserLoginView	user:signin	
    /user/update-affiliation	etalia.users.views.UserAffiliationUpdateView	user:update-affiliation	
    /user/update-email-digest-settings	etalia.users.views.UserEmailDigestSettingsUpdateView	user:update-email-digest-settings	
    /user/update-fingerprint-settings	etalia.users.views.UserFingerprintSettingsUpdateView	user:update-fingerprint-settings	
    /user/update-library	etalia.users.views.update_library	user:update-library	login_required
    /user/update-name	etalia.users.views.UpdateUserNameView	user:update-name	
    /user/update-position	etalia.users.views.UpdateUserPositionView	user:update-position	
    /user/update-stream-settings	etalia.users.views.UserStreamSettingsUpdateView	user:update-stream-settings	
    /user/update-title	etalia.users.views.UpdateUserTitleView	user:update-title	
    /user/update-trend-settings	etalia.users.views.UserTrendSettingsUpdateView	user:update-trend-settings	
    /user/user-init-check	etalia.users.views.user_init_check	user:user-init-check	login_required
    /user/user-update-settings-check	etalia.users.views.user_update_settings_check	user:user-update-settings-check	login_required
    /user/user-update-stream-check	etalia.users.views.user_update_stream_check	user:user-update-stream-check	login_required
    /user/user-update-trend-check	etalia.users.views.user_update_trend_check	user:user-update-trend-check	login_required